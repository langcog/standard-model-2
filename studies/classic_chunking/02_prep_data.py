"""Phonemize CHILDES caregiver utterances and CDI items with CMUdict (ARPAbet,
stress stripped, first pronunciation only), for CLASSIC input."""
import re
import sys
import pickle
from collections import Counter
import pandas as pd

W = sys.argv[1]            # work dir
CMU = sys.argv[2]          # cmudict.dict

pron = {}
with open(CMU, encoding='latin-1') as f:
    for line in f:
        parts = line.split('#')[0].strip().split()   # drop trailing "# name" etc. comments
        if not parts:
            continue
        w, ph = parts[0], parts[1:]
        if re.search(r'\(\d+\)$', w):
            continue                      # keep first pronunciation only
        pron[w] = tuple(re.sub(r'\d', '', p) for p in ph)
# hand-added pronunciations for frequent out-of-dictionary tokens (cf. the chapter's
# manual transcriptions for words with frequency >= 100)
MANUAL = {
    'uhhuh': 'AH HH AH', 'uhoh': 'AH OW', 'uhuh': 'AH AH', 'peekaboo': 'P IY K AH B UW',
    'horsie': 'HH AO R S IY', 'needta': 'N IY D T AH', 'sposta': 'S P OW S T AH',
    'useta': 'Y UW S T AH', 'hadta': 'HH AE D T AH', 'oink': 'OY N K', 'whee': 'W IY',
    'oo': 'UW', 'haha': 'HH AA HH AA', 'dyou': 'D Y UW', 'mkay': 'M EY K EY',
    'cockadoodledoo': 'K AA K AH D UW D AH L D UW', 'pattycake': 'P AE T IY K EY K',
    'grrr': 'G ER R', 'mmhm': 'M HH M', 'hmm': 'HH M', 'mm': 'M', 'mmm': 'M',
}
for k, v in MANUAL.items():
    pron.setdefault(k, tuple(v.split()))
phonemes = sorted({p for v in pron.values() for p in v})
print('cmudict entries', len(pron), 'phonemes', len(phonemes))


def split_compound(t):
    """OOV token -> two dictionary words (e.g. icecream, highchair), else None."""
    for i in range(3, len(t) - 2):
        if t[:i] in pron and t[i:] in pron:
            return [t[:i], t[i:]]
    return None


def tokens(text):
    text = text.lower().replace('_', ' ').replace('+', ' ').replace('-', ' ')
    text = re.sub(r"[^a-z' ]", ' ', text)
    out = []
    for t in text.split():
        t = t.strip("'") if t.strip("'") in pron else t
        if not t:
            continue
        if t not in pron:
            sp = split_compound(t)
            if sp:
                out.extend(sp)
                continue
        out.append(t)
    return out


# ---- CDI items -> phoneme-string variants -----------------------------------
items = pd.read_csv(f'{W}/v2_items_delta.csv')
SKIP = {"babysitter's name", "child's own name", "pet's name"}
targets, dropped = {}, []
for r in items.itertuples():
    d = r.item_definition
    if d in SKIP:
        dropped.append(d)
        continue
    variants = []
    for v in d.split('/'):
        v = re.sub(r'\(.*?\)', '', v).replace('*', '').replace('!', '').strip().lower()
        ws = tokens(v)
        if ws and all(w in pron for w in ws):
            variants.append(' '.join(' '.join(pron[w]) for w in ws))
    if variants:
        targets[int(r.jj)] = sorted(set(variants))
    else:
        dropped.append(d)
print('CDI items with phonemes:', len(targets), '/', len(items), 'dropped:', dropped)

# ---- CHILDES utterances -----------------------------------------------------
u = pd.read_csv(f'{W}/childes_eng_na_utts.csv.gz')
print('utterances read', len(u))
oov = Counter()
utts, ages, keep = [], [], []
for g, age in zip(u.gloss.astype(str), u.target_child_age):
    ws = tokens(g)
    if not ws or any(w in ('xxx', 'yyy', 'www') for w in ws):
        keep.append(False)
        continue
    bad = [w for w in ws if w not in pron]
    if bad:
        oov.update(bad)
        keep.append(False)
        continue
    utts.append(tuple(pron[w] for w in ws))
    ages.append(age)
    keep.append(True)
print('kept', len(utts), f'({len(utts) / len(u):.1%}); top OOV:', oov.most_common(40))
with open(f'{W}/prep.pkl', 'wb') as f:
    pickle.dump({'utts': utts, 'ages': ages, 'targets': targets, 'phonemes': phonemes,
                 'oov': dict(oov), 'dropped_items': dropped}, f)
