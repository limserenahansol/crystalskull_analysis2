# Push to GitHub

To upload this code to https://github.com/limserenahansol/crystalskull_analysis:

## Option A: Replace existing repo content (recommended)

The crystalskull_analysis repo currently has behavior_task code. To replace it with this crystal skull pipeline:

```powershell
cd C:\Users\hsollim\Desktop\crystalskull_analysis_repo

# Add remote
git remote add origin https://github.com/limserenahansol/crystalskull_analysis.git

# Rename branch to main (GitHub default)
git branch -M main

# Force push to replace remote content
git push -u origin main --force
```

**Warning:** `--force` overwrites the remote. Use only if you intend to replace the current repo content.

## Option B: Push to a new branch

If you want to keep the existing content and add this as a branch:

```powershell
cd C:\Users\hsollim\Desktop\crystalskull_analysis_repo

git remote add origin https://github.com/limserenahansol/crystalskull_analysis.git
git branch -M crystal_skull_analysis
git push -u origin crystal_skull_analysis
```

Then merge or switch branches on GitHub as needed.

---

After pushing, you can delete this `PUSH_INSTRUCTIONS.md` file or remove it from the repo.
