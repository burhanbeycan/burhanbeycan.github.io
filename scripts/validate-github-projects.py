#!/usr/bin/env python3
"""Validate the GitHub portfolio project data file.

This script is intentionally dependency-free so it can run locally and in GitHub Actions.
It checks the structure of assets/data/github-projects.json and fails fast when required
fields are missing, duplicated, or inconsistent.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

PROJECT_REQUIRED_FIELDS = {
    "id",
    "title",
    "category",
    "status",
    "priority",
    "repository",
    "github_url",
    "live_url",
    "tagline",
    "why_it_matters",
    "technical_stack",
    "deliverables",
    "evidence_to_show",
}

ALLOWED_STATUSES = {"active_repo", "planned", "planned_or_private_work"}
SLUG_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
URL_RE = re.compile(r"^https://")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing data file: {path}")

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def validate_profile(data: dict[str, Any], errors: list[str]) -> None:
    profile = data.get("profile")
    require(isinstance(profile, dict), "profile must be an object", errors)
    if not isinstance(profile, dict):
        return

    for field in ["name", "positioning", "website", "github_owner"]:
        require(bool(profile.get(field)), f"profile.{field} is required", errors)

    if profile.get("website"):
        require(bool(URL_RE.match(profile["website"])), "profile.website must start with https://", errors)


def validate_categories(data: dict[str, Any], errors: list[str]) -> set[str]:
    categories = data.get("categories")
    require(isinstance(categories, list) and categories, "categories must be a non-empty list", errors)
    category_ids: set[str] = set()

    if not isinstance(categories, list):
        return category_ids

    for index, category in enumerate(categories):
        require(isinstance(category, dict), f"categories[{index}] must be an object", errors)
        if not isinstance(category, dict):
            continue

        category_id = category.get("id")
        require(bool(category_id), f"categories[{index}].id is required", errors)
        if category_id:
            require(category_id not in category_ids, f"duplicate category id: {category_id}", errors)
            require(bool(SLUG_RE.match(category_id)), f"category id must be a slug: {category_id}", errors)
            category_ids.add(category_id)

        for field in ["label", "summary"]:
            require(bool(category.get(field)), f"categories[{index}].{field} is required", errors)

    return category_ids


def validate_projects(data: dict[str, Any], category_ids: set[str], errors: list[str]) -> None:
    projects = data.get("projects")
    require(isinstance(projects, list) and projects, "projects must be a non-empty list", errors)
    if not isinstance(projects, list):
        return

    project_ids: set[str] = set()
    priorities: set[int] = set()

    for index, project in enumerate(projects):
        require(isinstance(project, dict), f"projects[{index}] must be an object", errors)
        if not isinstance(project, dict):
            continue

        missing = sorted(PROJECT_REQUIRED_FIELDS - project.keys())
        require(not missing, f"projects[{index}] is missing fields: {', '.join(missing)}", errors)

        project_id = project.get("id")
        if project_id:
            require(project_id not in project_ids, f"duplicate project id: {project_id}", errors)
            require(bool(SLUG_RE.match(project_id)), f"project id must be a slug: {project_id}", errors)
            project_ids.add(project_id)

        category = project.get("category")
        require(category in category_ids, f"project {project_id or index} has unknown category: {category}", errors)

        status = project.get("status")
        require(status in ALLOWED_STATUSES, f"project {project_id or index} has invalid status: {status}", errors)

        priority = project.get("priority")
        require(isinstance(priority, int), f"project {project_id or index} priority must be an integer", errors)
        if isinstance(priority, int):
            require(priority > 0, f"project {project_id or index} priority must be positive", errors)
            require(priority not in priorities, f"duplicate project priority: {priority}", errors)
            priorities.add(priority)

        for field in ["title", "repository", "tagline", "why_it_matters"]:
            require(bool(project.get(field)), f"project {project_id or index}.{field} is required", errors)

        for url_field in ["github_url", "live_url"]:
            url = project.get(url_field, "")
            if url:
                require(bool(URL_RE.match(url)), f"project {project_id or index}.{url_field} must start with https://", errors)

        for list_field in ["technical_stack", "deliverables", "evidence_to_show"]:
            value = project.get(list_field)
            require(isinstance(value, list) and value, f"project {project_id or index}.{list_field} must be a non-empty list", errors)
            if isinstance(value, list):
                require(all(isinstance(item, str) and item.strip() for item in value), f"project {project_id or index}.{list_field} must contain non-empty strings", errors)


def validate(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        data = load_json(path)
    except Exception as exc:  # noqa: BLE001 - command-line validator should report all load errors simply.
        return [str(exc)]

    validate_profile(data, errors)
    category_ids = validate_categories(data, errors)
    validate_projects(data, category_ids, errors)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate GitHub portfolio project data.")
    parser.add_argument(
        "data_file",
        nargs="?",
        default="assets/data/github-projects.json",
        help="Path to the JSON data file. Default: assets/data/github-projects.json",
    )
    args = parser.parse_args()

    errors = validate(Path(args.data_file))
    if errors:
        print("GitHub project data validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print(f"GitHub project data validation passed: {args.data_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
