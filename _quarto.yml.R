project:
  type: book

book:
  title: "サイクリングとウォーキングの健康効果比較"
author: "Kizen Sasaki"
date-format: 'YYYY MMMM DD'
date: 2024-07-20
date-modified: today
chapters:
- index.qmd
- 01-intro.qmd
- 02-methods.qmd
- 03-results.qmd
- 04-discussion.qmd
- PRISMA.qmd
- references.qmd

bibliography: references.bib

format:
  html:
  number-sections: true
number-depth: 1
theme: cosmo
pdf:
  documentclass: scrreprt

crossref:
  fig-prefix: "図"
fig-title: "図"
fig-labels: arabic

language:
  label:
  fig: "図 "
tab: "表 "