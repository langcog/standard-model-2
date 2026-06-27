"""Grammar eval for CHILDES GPT-2: BLiMP (67 paradigms, HF) + Zorro (23 paradigms, CHILDES-vocab).

Minimal-pair scoring: for each (good, bad) pair, the model is "correct" if it assigns
higher total log-probability (lower cross-entropy) to the grammatical sentence. This is the
standard BLiMP/Zorro autoregressive-LM scoring. Sentences in a pair are length-matched, so
we use summed token log-prob (no length normalization), matching the convention.

Zorro convention (from the repo README/data.py): in each paradigm file, 0-based even lines
(1-indexed odd) are UNgrammatical, 0-based odd lines (1-indexed even) are grammatical.

Usable two ways:
  - standalone:  python eval_grammar.py --model_dir D --tokenizer_dir T --zorro_dir Z --out out.json
  - imported:    from eval_grammar import grammar_accuracy; grammar_accuracy(model, tok, ...)
"""
import argparse, glob, json, os
import torch
import torch.nn.functional as F


def _batch_sent_logprob(model, tokenizer, sents, device, batch_size, max_len=128):
    """Sum of per-token log-probs for each sentence (autoregressive). Returns list[float]."""
    out = []
    model.eval()
    for i in range(0, len(sents), batch_size):
        chunk = sents[i:i + batch_size]
        enc = tokenizer(chunk, return_tensors="pt", padding=True, truncation=True, max_length=max_len)
        ids = enc["input_ids"].to(device)
        attn = enc["attention_mask"].to(device)
        with torch.no_grad():
            logits = model(input_ids=ids, attention_mask=attn).logits  # [B,T,V]
            logp = F.log_softmax(logits.float(), dim=-1)
            tgt = ids[:, 1:]                                  # predict token t from <t
            pred = logp[:, :-1, :]
            tok_lp = pred.gather(-1, tgt.unsqueeze(-1)).squeeze(-1)   # [B,T-1]
            m = attn[:, 1:].float()
            sent_lp = (tok_lp * m).sum(dim=1)                 # [B]
        out.extend(sent_lp.tolist())
    return out


def _accuracy_from_pairs(model, tokenizer, good, bad, device, batch_size, max_len):
    """good[i] vs bad[i] minimal pairs -> accuracy (model prefers good)."""
    lp_good = _batch_sent_logprob(model, tokenizer, good, device, batch_size, max_len)
    lp_bad = _batch_sent_logprob(model, tokenizer, bad, device, batch_size, max_len)
    correct = sum(1 for g, b in zip(lp_good, lp_bad) if g > b)
    return correct / max(1, len(good)), len(good)


def eval_zorro(model, tokenizer, zorro_dir, device, batch_size=64, max_pairs=0, max_len=128):
    files = sorted(glob.glob(os.path.join(zorro_dir, "sentences", "babyberta", "*.txt")))
    per = {}
    for f in files:
        para = os.path.splitext(os.path.basename(f))[0]
        lines = [ln.strip() for ln in open(f) if ln.strip()]
        bad, good = lines[0::2], lines[1::2]          # 0-based even=bad, odd=good
        n = min(len(good), len(bad))
        good, bad = good[:n], bad[:n]
        if max_pairs and n > max_pairs:
            good, bad = good[:max_pairs], bad[:max_pairs]
        acc, npair = _accuracy_from_pairs(model, tokenizer, good, bad, device, batch_size, max_len)
        per[para] = {"acc": acc, "n": npair}
    overall = sum(v["acc"] for v in per.values()) / max(1, len(per))   # paradigm-macro avg
    return {"overall": overall, "n_paradigms": len(per), "per_paradigm": per}


# The 67 BLiMP paradigms (stable set). Hardcoded so we don't call the HF registry
# (get_dataset_config_names) -- under HF_DATASETS_OFFLINE it returns 'default' and breaks,
# and online it rate-limits when many arms hit it at once.
BLIMP_CONFIGS = [
    "adjunct_island", "anaphor_gender_agreement", "anaphor_number_agreement",
    "animate_subject_passive", "animate_subject_trans", "causative", "complex_NP_island",
    "coordinate_structure_constraint_complex_left_branch",
    "coordinate_structure_constraint_object_extraction", "determiner_noun_agreement_1",
    "determiner_noun_agreement_2", "determiner_noun_agreement_irregular_1",
    "determiner_noun_agreement_irregular_2", "determiner_noun_agreement_with_adj_2",
    "determiner_noun_agreement_with_adj_irregular_1", "determiner_noun_agreement_with_adj_irregular_2",
    "determiner_noun_agreement_with_adjective_1", "distractor_agreement_relational_noun",
    "distractor_agreement_relative_clause", "drop_argument", "ellipsis_n_bar_1", "ellipsis_n_bar_2",
    "existential_there_object_raising", "existential_there_quantifiers_1",
    "existential_there_quantifiers_2", "existential_there_subject_raising",
    "expletive_it_object_raising", "inchoative", "intransitive",
    "irregular_past_participle_adjectives", "irregular_past_participle_verbs",
    "irregular_plural_subject_verb_agreement_1", "irregular_plural_subject_verb_agreement_2",
    "left_branch_island_echo_question", "left_branch_island_simple_question",
    "matrix_question_npi_licensor_present", "npi_present_1", "npi_present_2",
    "only_npi_licensor_present", "only_npi_scope", "passive_1", "passive_2",
    "principle_A_c_command", "principle_A_case_1", "principle_A_case_2", "principle_A_domain_1",
    "principle_A_domain_2", "principle_A_domain_3", "principle_A_reconstruction",
    "regular_plural_subject_verb_agreement_1", "regular_plural_subject_verb_agreement_2",
    "sentential_negation_npi_licensor_present", "sentential_negation_npi_scope",
    "sentential_subject_island", "superlative_quantifiers_1", "superlative_quantifiers_2",
    "tough_vs_raising_1", "tough_vs_raising_2", "transitive", "wh_island",
    "wh_questions_object_gap", "wh_questions_subject_gap", "wh_questions_subject_gap_long_distance",
    "wh_vs_that_no_gap", "wh_vs_that_no_gap_long_distance", "wh_vs_that_with_gap",
    "wh_vs_that_with_gap_long_distance",
]


def eval_blimp(model, tokenizer, device, batch_size=64, max_pairs=0, max_len=128, configs=None):
    from datasets import load_dataset
    configs = configs or BLIMP_CONFIGS
    per = {}
    for cfg in configs:
        try:
            ds = load_dataset("blimp", cfg, split="train")
        except Exception as e:
            print(f"[grammar] skip blimp/{cfg}: {repr(e)[:90]}", flush=True)
            continue
        good = list(ds["sentence_good"]); bad = list(ds["sentence_bad"])
        if max_pairs and len(good) > max_pairs:
            good, bad = good[:max_pairs], bad[:max_pairs]
        acc, npair = _accuracy_from_pairs(model, tokenizer, good, bad, device, batch_size, max_len)
        per[cfg] = {"acc": acc, "n": npair}
    overall = sum(v["acc"] for v in per.values()) / max(1, len(per))
    return {"overall": overall, "n_paradigms": len(per), "per_paradigm": per}


def grammar_accuracy(model, tokenizer, zorro_dir, device, batch_size=64,
                     max_pairs=0, max_len=128, do_blimp=True, blimp_configs=None):
    """Convenience wrapper used by the training callback and CLI. Returns dict with
    zorro/blimp overall accuracies (paradigm-macro-averaged)."""
    res = {"zorro": eval_zorro(model, tokenizer, zorro_dir, device, batch_size, max_pairs, max_len)}
    if do_blimp:
        res["blimp"] = eval_blimp(model, tokenizer, device, batch_size, max_pairs, max_len, blimp_configs)
    res["zorro_overall"] = res["zorro"]["overall"]
    res["blimp_overall"] = res.get("blimp", {}).get("overall")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--zorro_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--max_pairs", type=int, default=0, help="subsample pairs/paradigm (0=all)")
    ap.add_argument("--max_len", type=int, default=128)
    ap.add_argument("--no_blimp", action="store_true")
    args = ap.parse_args()

    from transformers import AutoTokenizer, GPT2LMHeadModel
    device = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(args.tokenizer_dir)
    if tok.pad_token is None:
        tok.pad_token = tok.eos_token
    model = GPT2LMHeadModel.from_pretrained(args.model_dir).to(device)
    res = grammar_accuracy(model, tok, args.zorro_dir, device, args.batch_size,
                           args.max_pairs, args.max_len, do_blimp=not args.no_blimp)
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[grammar] zorro={res['zorro_overall']:.4f} "
          f"blimp={res['blimp_overall'] if res['blimp_overall'] is not None else 'NA'}  -> {args.out}")


if __name__ == "__main__":
    main()
