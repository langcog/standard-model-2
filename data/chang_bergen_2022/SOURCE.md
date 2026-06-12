# Chang & Bergen (2022) per-word sigmoid fits — provenance

The four `*_sigmoids.txt` files in this directory are **byte-for-byte copies**
of the per-word sigmoid fits released by Chang & Bergen (2022), *Word
Acquisition in Neural Language Models* (TACL).

## Upstream source

Repo: <https://github.com/tylerachang/word-acquisition-language-models>
Path: `r_code/tacl_data/lm_data/processed/{bert,gpt2,lstm,bilstm}_sigmoids.txt`
Raw:  `https://raw.githubusercontent.com/tylerachang/word-acquisition-language-models/main/r_code/tacl_data/lm_data/processed/<model>_sigmoids.txt`

MD5 (verified identical to upstream `main`, 2026-06):

| file                | md5                                |
|---------------------|------------------------------------|
| `bert_sigmoids.txt`   | `c2cb9829b4ba3e5f72407821cfc0ca33` |
| `bilstm_sigmoids.txt` | `778b7869a56a5baa3feaa8349f7c9fc7` |
| `gpt2_sigmoids.txt`   | `ab2fef4c2ba3b15a5a8621ff1966e701` |
| `lstm_sigmoids.txt`   | `bb2dc7ac5da1a98cabfa604fecd3beae` |

Re-download / re-verify:

```sh
base=https://raw.githubusercontent.com/tylerachang/word-acquisition-language-models/main/r_code/tacl_data/lm_data/processed
for m in bert bilstm gpt2 lstm; do curl -sO "$base/${m}_sigmoids.txt"; done
```

## What the files contain

One row per CDI word (611 of 651 words that are single tokens), per
architecture. Each row is a 4-parameter logistic (4-PL) fit to the model's
mean token surprisal trajectory over training steps (sampled at ~200 log-spaced
checkpoints). Columns:

`Token  MaxSurprisal  MinSurprisal  ParamUpper  ParamLower  ParamXmid  ParamScale`

The fitted curve is
`surprisal(x) = ParamLower + (ParamUpper - ParamLower) / (1 + exp((x - ParamXmid)/ParamScale))`,
with `x = log10(training steps)`. These are C&B's own fits (produced by their
`r_code/tacl_analyses.rmd`); we do **not** re-fit them.

## The four models (C&B §3.1, Table 1)

All four are single models (no seed replication), trained on the **same**
corpus: a combination of **BookCorpus (Zhu et al. 2015) + WikiText-103 (Merity
et al. 2017)**, 25.6M sentence pairs, SentencePiece unigram tokenizer, 1M steps.

| Architecture | # Params | Eval perplexity | Objective                       |
|--------------|----------|-----------------|---------------------------------|
| LSTM         | 37M      | 54.8            | unidirectional, 3 stacked layers|
| BiLSTM       | 51M      | 9.0             | bidirectional LSTM              |
| GPT-2        | 108M     | 30.2            | causal LM (GPT-2 size)          |
| BERT         | 109M     | 7.2             | masked LM (BERT-base)           |

## How we use them

`paper/build_cache.R` §6 converts each per-word `ParamScale` to a slope on
**natural-log** experience, comparable to children's kappa_i:

```
slope_natural = (1 / ParamScale) / ln(10)
```

Degenerate fits are dropped (`ParamScale` in (0.01, 10) and surprisal range
`ParamUpper - ParamLower > 1`). All four architectures are pooled into one
"LMs: C&B 2022 (4 architectures)" population (~2,410 word fits) for the
acceleration figure (`fig-llm-acceleration`). Pooled median slope ≈ 0.84 —
i.e. these models sit at the kappa = 1 unit-accumulator boundary, ~10x below
children.
