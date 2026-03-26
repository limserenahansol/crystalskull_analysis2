# Pushing to GitHub

**Canonical repo:** [limserenahansol/crystalskull_analysis2](https://github.com/limserenahansol/crystalskull_analysis2)

From this folder:

```bash
git remote -v
git add -A
git commit -m "Your message"
git push origin main
```

If `origin` still points at an old or renamed repository:

```bash
git remote set-url origin https://github.com/limserenahansol/crystalskull_analysis2.git
git push origin main
```
