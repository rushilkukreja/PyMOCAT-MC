# Test Paper Compilation Locally

## Option 1: Use Docker (Most Reliable)
```bash
# Pull the JOSS Docker image
docker pull openjournals/paperdraft

# Run from your repository root
docker run --rm \
    --volume $PWD/paper:/app \
    --env JOURNAL=joss \
    openjournals/paperdraft
```

## Option 2: Install Dependencies Locally
```bash
# Install pandoc and dependencies
brew install pandoc pandoc-citeproc  # macOS
# or
sudo apt-get install pandoc pandoc-citeproc  # Ubuntu

# Try basic compilation
cd paper
pandoc paper.md --bibliography=paper.bib --citeproc -o test-output.pdf
```

## Option 3: Minimal GitHub Action Test
Create a simpler action that just validates the paper structure:

```yaml
name: Validate Paper
on: workflow_dispatch
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: |
          echo "Paper exists: $(test -f paper/paper.md && echo 'YES' || echo 'NO')"
          echo "Bib exists: $(test -f paper/paper.bib && echo 'YES' || echo 'NO')"
          head -30 paper/paper.md
```

## What to try first:
1. **Docker approach** (if you have Docker installed)
2. **Simple validation action** (just to test GitHub Actions work)
3. **Full PDF generation action** (after confirming the basics work)