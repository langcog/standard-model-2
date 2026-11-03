"""
Fast re-implementation of CLASSIC (Chunking Lexical and Sub-lexical Sequences in
Children; Jones et al.), following the OSF script `classic.py` of Rowland et al.
(2025, Language Learning; https://osf.io/9s4d7) rule-for-rule, but with set/dict
lookups instead of pandas scans.

Faithful to the OSF script:
  * chunk access: chunk with recency x (x=0 is the utterance-final chunk) is
    accessed with prob 1 / (1 + exp(-(L - x - 1))), L = 4.5
  * learning: each adjacent pair of accessed chunks is learned with prob LR = 0.10
  * parse order per utterance: multiword -> lexical (whole word) -> sublexical
    (greedy longest, left to right, within a word) -> single phonemes
  * pair learning: same word -> lexical if that word has exactly two chunks,
    else sublexical; across words -> multiword iff both are lexical/multiword
  * second pass: an accessed chunk (learning_ok again) is also stored under the
    type it was parsed as, if not yet stored under that type
  * a vocabulary item is "known" iff its phoneme string is a lexical or
    multiword chunk

Differences (deliberate): lexicon is a set of strings (the OSF script appends
duplicate rows); counters (n_use, n_access, ...) that never affect behaviour are
dropped.
"""
import math
import random

ACCESS_L = 4.5
LEARNING_RATE = 0.10


def _access_table(L, n=60):
    return [1.0 / (1.0 + math.exp(-(L - x - 1))) for x in range(n)]


def parse(words, allch, maxlen, acc_p, rnd):
    """words: list of tuples of phoneme strings. Returns list of chunk tuples
    (string, type, word_index, accessed)."""
    nW = len(words)
    flat = [p for w in words for p in w]
    nP = len(flat)
    wstart = []
    c = 0
    for w in words:
        wstart.append(c)
        c += len(w)
    cid_of = [-1] * nP
    ctype = [None] * nP
    cid = 0

    # multiword
    s = 0
    e = nW - 1
    while e - s > 0:
        a = wstart[s]
        b = wstart[e] + len(words[e])
        if b - a <= maxlen and ' '.join(flat[a:b]) in allch:
            for i in range(a, b):
                cid_of[i] = cid
                ctype[i] = 'multiword'
            cid += 1
            s = e + 1
            e = nW - 1
        else:
            e -= 1
            if e == s:
                s += 1
                e = nW - 1

    # lexical
    for wp in range(nW):
        a = wstart[wp]
        b = a + len(words[wp])
        if cid_of[a] == -1 and ' '.join(flat[a:b]) in allch:
            for i in range(a, b):
                cid_of[i] = cid
                ctype[i] = 'lexical'
            cid += 1

    # sublexical (greedy longest-first within an unassigned word)
    for wp in range(nW):
        wa = wstart[wp]
        wb = wa + len(words[wp])
        if cid_of[wa] != -1:
            continue
        s = wa
        e = wb
        while e - s > 1:
            if ' '.join(flat[s:e]) in allch:
                for i in range(s, e):
                    cid_of[i] = cid
                    ctype[i] = 'sublexical'
                cid += 1
                s = e
                e = wb
            else:
                e -= 1
                if e == s + 1:
                    s += 1
                    e = wb

    # phonemic
    for i in range(nP):
        if cid_of[i] == -1:
            cid_of[i] = cid
            ctype[i] = 'phonemic'
            cid += 1

    # chunks in positional order
    word_of = []
    for wp in range(nW):
        word_of.extend([wp] * len(words[wp]))
    chunks = []
    i = 0
    while i < nP:
        j = i
        while j + 1 < nP and cid_of[j + 1] == cid_of[i]:
            j += 1
        chunks.append([' '.join(flat[i:j + 1]), ctype[i], word_of[i], False])
        i = j + 1

    n = len(chunks)
    for k in range(n):
        x = n - k - 1
        p = acc_p[x] if x < len(acc_p) else 0.0
        chunks[k][3] = rnd() <= p
    return chunks


def learn(chunks, lex, multi, sub, allch, rnd, lr, on_word):
    """Update lexicon in place; on_word(string) is called when a string first
    enters lex or multi."""
    n = len(chunks)
    per_word = {}
    for ch in chunks:
        per_word[ch[2]] = per_word.get(ch[2], 0) + 1
    new = []  # (string, type)
    for i in range(n - 1):
        a, b = chunks[i], chunks[i + 1]
        if a[3] and b[3] and rnd() <= lr:
            s = a[0] + ' ' + b[0]
            if a[2] == b[2]:
                new.append((s, 'lexical' if per_word[a[2]] == 2 else 'sublexical'))
            elif a[1] in ('lexical', 'multiword') and b[1] in ('lexical', 'multiword'):
                new.append((s, 'multiword'))
    for ch in chunks:
        if ch[3] and rnd() <= lr:
            t = ch[1]
            if t == 'sublexical' and ch[0] not in sub:
                new.append((ch[0], 'sublexical'))
            elif t == 'lexical' and ch[0] not in lex:
                new.append((ch[0], 'lexical'))
            elif t == 'multiword' and ch[0] not in multi:
                new.append((ch[0], 'multiword'))
    grew_max = 0
    for s, t in new:
        allch.add(s)
        grew_max = max(grew_max, s.count(' ') + 1)
        if t == 'lexical':
            if s not in lex:
                lex.add(s)
                if s not in multi:
                    on_word(s)
        elif t == 'multiword':
            if s not in multi:
                multi.add(s)
                if s not in lex:
                    on_word(s)
        else:
            sub.add(s)
    return grew_max


def run(utts, targets, seed, phonemes, lr=LEARNING_RATE, L=ACCESS_L,
        n_checkpoints=20):
    """utts: list of utterances (list of word tuples of phonemes), in input order.
    targets: dict item_id -> list of phoneme-string variants (each a space-joined
    string). Returns dict with per-item learn index and exposure counts, plus
    vocabulary growth at checkpoints."""
    rng = random.Random(seed)
    rnd = rng.random
    acc_p = _access_table(L)
    allch = set(phonemes)
    lex, multi, sub = set(), set(), set()
    maxlen = 1
    str2items = {}
    for item, variants in targets.items():
        for v in variants:
            str2items.setdefault(v, []).append(item)
    learned_at = {}
    n = len(utts)
    cps = {int(round(n * (k + 1) / n_checkpoints)) for k in range(n_checkpoints)}
    growth = []
    cur = [0]

    def on_word(s):
        for item in str2items.get(s, ()):
            if item not in learned_at:
                learned_at[item] = cur[0]

    for u, words in enumerate(utts, start=1):
        cur[0] = u
        chunks = parse(words, allch, maxlen, acc_p, rnd)
        m = learn(chunks, lex, multi, sub, allch, rnd, lr, on_word)
        if m > maxlen:
            maxlen = m  # longest stored chunk, in phonemes (prunes multiword search)
        if u in cps:
            known = sum(1 for item in targets if item in learned_at)
            growth.append((u, known, len(lex), len(multi), len(sub)))
    return {'learned_at': learned_at, 'growth': growth, 'n_utts': n,
            'sizes': (len(lex), len(multi), len(sub))}
