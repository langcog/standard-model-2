#!/usr/bin/env python
"""Export a trained GPT-2 checkpoint -> bf16 safetensors, upload to HF (+ optional Oak).
Streaming archival: caller deletes the local checkpoint after this returns OK."""
import argparse, os, shutil, subprocess
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer
from huggingface_hub import HfApi

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)          # output_dir/checkpoint-N (best, save_total_limit=1)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--repo", required=True)           # mcxfrank/childes-gpt2-ladder
    ap.add_argument("--path_in_repo", required=True)   # e.g. development/seed0/rung0.5M
    ap.add_argument("--export_dir", required=True)     # small scratch (~250MB), removed at end
    ap.add_argument("--meta", default="")              # provenance JSON string
    ap.add_argument("--oak", default="")               # optional oak-dtn:.../path
    a = ap.parse_args()

    if os.path.exists(a.export_dir):
        shutil.rmtree(a.export_dir)
    os.makedirs(a.export_dir)

    # best checkpoint -> bf16 (training compute dtype) -> self-contained model dir
    model = AutoModelForCausalLM.from_pretrained(a.ckpt, torch_dtype=torch.bfloat16)
    model.save_pretrained(a.export_dir, safe_serialization=True)
    AutoTokenizer.from_pretrained(a.tokenizer_dir).save_pretrained(a.export_dir)
    if a.meta:
        with open(os.path.join(a.export_dir, "provenance.json"), "w") as f:
            f.write(a.meta)

    api = HfApi()
    api.upload_folder(folder_path=a.export_dir, repo_id=a.repo,
                      path_in_repo=a.path_in_repo, repo_type="model",
                      commit_message="add " + a.path_in_repo)

    if a.oak:
        subprocess.run(["rsync", "-rlt", "--partial", a.export_dir + "/", a.oak], check=True)

    shutil.rmtree(a.export_dir)
    print("EXPORT_UPLOAD_OK", a.path_in_repo, flush=True)

if __name__ == "__main__":
    main()
