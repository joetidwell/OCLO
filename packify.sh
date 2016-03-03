#!/bin/bash

echo "Delete old documentation files"
rm NAMESPACE
rm man/*

echo "Create new documentation"
Rscript -e "library(devtools); library(roxygen2);document()"

echo "Commit changes to git repo"
git add *
git commit -m "Repackage and document"

echo "Push to gitHub repo"
git push
