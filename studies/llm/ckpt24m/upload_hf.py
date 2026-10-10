"""Upload the gated 24M trajectories to mcxfrank/childes-gpt2-ladder under trajectory/,
then insert the trajectory section into the model card. Run on ccn2 from ckpt24m/.
Staging uses hard links (same filesystem, no extra disk)."""
import os
import shutil
import subprocess

from huggingface_hub import HfApi, hf_hub_download

REPO, H = "mcxfrank/childes-gpt2-ladder", os.path.dirname(os.path.abspath(__file__))
stage = os.path.join(H, "hf_stage")
shutil.rmtree(stage, ignore_errors=True)
for s in (42, 0, 123):
    src, dst = f"{H}/runs/seed{s}/trajectory", f"{stage}/trajectory/seed{s}/rung24M"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    subprocess.run(["cp", "-al", src, dst], check=True)
    os.link(f"{H}/runs/seed{s}/traj_surprisal.csv", f"{dst}/surprisal.csv")

api = HfApi()
api.upload_large_folder(repo_id=REPO, repo_type="model", folder_path=stage)

card = open(hf_hub_download(REPO, "README.md")).read()
section = open(f"{H}/README_trajectory_section.md").read()
if "## Training trajectories" not in card:
    card = card.replace("## Load a model", section + "## Load a model", 1)
    api.upload_file(path_or_fileobj=card.encode(), path_in_repo="README.md", repo_id=REPO,
                    commit_message="model card: document trajectory/ checkpoints")
n = sum(f.startswith("trajectory/") for f in api.list_repo_files(REPO))
print("UPLOAD_OK", n, "files under trajectory/")
