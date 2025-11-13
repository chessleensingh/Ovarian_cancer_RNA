#!/bin/bash

# Quick push script for Ovarian Cancer RNA analysis
# Run this script to push to GitHub

cd "$(dirname "$0")"

echo "==========================================="
echo "Pushing to GitHub: Ovarian_cancer_RNA"
echo "==========================================="
echo ""

# Show what will be pushed
echo "Files to be pushed:"
git status --short
echo ""

# Push to GitHub
echo "Pushing to GitHub..."
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Successfully pushed to GitHub!"
    echo ""
    echo "View your repository at:"
    echo "https://github.com/chessleensingh/Ovarian_cancer_RNA"
    echo ""
    echo "⚠️  IMPORTANT: Revoke your GitHub token now!"
    echo "https://github.com/settings/tokens"
else
    echo ""
    echo "✗ Push failed. Try running manually:"
    echo "cd ~/Desktop/Ovarian_Cancer_GitHub"
    echo "git push -u origin main"
fi
