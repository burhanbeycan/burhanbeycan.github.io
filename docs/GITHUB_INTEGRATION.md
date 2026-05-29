# GitHub Portfolio Integration Guide

This guide explains how the GitHub-integrated project portfolio is structured and how to maintain it.

## Files added

- `github-projects.html`: standalone project page.
- `assets/data/github-projects.json`: single source of truth for project cards.
- `assets/js/github-projects.js`: renders project cards, filters, category summaries, and GitHub metadata.
- `assets/css/github-projects.css`: page-specific styling.
- `scripts/validate-github-projects.py`: local validation script for the JSON data.
- `.github/workflows/validate-github-projects.yml`: automated validation workflow.
- `docs/PROJECT_PORTFOLIO_STRATEGY.md`: positioning and execution roadmap.

## Editing project data

Edit this file when you want to add or change a project:

```text
assets/data/github-projects.json
```

Each project should include an ID, title, category, status, priority, repository name, GitHub URL, optional live URL, tagline, career signal, technical stack, deliverables, and evidence to show.

## Current categories

- `general-ai`
- `computer-vision`
- `scientific-ai`
- `research-tools`

Keep the categories focused. The purpose is to show a balanced AI portfolio rather than a long, unfocused project list.

## Current statuses

- `active_repo`
- `planned`
- `planned_or_private_work`

Use `active_repo` only when a repository is public and useful enough to show to a reviewer.

## Validation

Run this before committing changes:

```bash
python scripts/validate-github-projects.py assets/data/github-projects.json
```

The same check runs automatically through GitHub Actions for relevant pull requests and pushes to `main`.

## Repository quality checklist

Before a repository becomes a strong portfolio item, add:

- A concise README explaining the problem, method, and result.
- Setup instructions.
- One screenshot, diagram, or demo video.
- Data and model limitations.
- Evaluation metrics or a baseline comparison.
- A short paragraph explaining the career signal.

## Suggested homepage link

Add a navigation link from the main homepage to `github-projects.html`. Recommended placement is after the existing Research item and before Publications.

Also add a GitHub profile icon near the other social links so reviewers can reach the GitHub profile quickly.
