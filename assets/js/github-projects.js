(() => {
  const DATA_URL = 'assets/data/github-projects.json';
  const OWNER_FALLBACK = 'burhanbeycan';
  const state = {
    data: null,
    repoMeta: new Map(),
    activeFilter: 'all'
  };

  const els = {
    filterRow: document.querySelector('[data-project-filters]'),
    projectGrid: document.querySelector('[data-project-grid]'),
    summaryGrid: document.querySelector('[data-summary-grid]'),
    updated: document.querySelector('[data-last-updated]')
  };

  const escapeHtml = (value = '') => String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#039;');

  const statusLabel = (status = '') => {
    const labels = {
      active_repo: 'Active repo',
      planned: 'Planned build',
      planned_or_private_work: 'Planned/private'
    };
    return labels[status] || status.replace(/_/g, ' ');
  };

  const formatDate = (value) => {
    if (!value) return '';
    try {
      return new Intl.DateTimeFormat('en', { month: 'short', year: 'numeric' }).format(new Date(value));
    } catch (_) {
      return '';
    }
  };

  const categoryById = () => new Map((state.data?.categories || []).map((category) => [category.id, category]));

  async function getJson(url) {
    const response = await fetch(url, { cache: 'no-store' });
    if (!response.ok) throw new Error(`Could not load ${url}: ${response.status}`);
    return response.json();
  }

  async function getRepoMeta(owner, repo) {
    if (!repo || state.repoMeta.has(repo)) return state.repoMeta.get(repo) || null;

    try {
      const response = await fetch(`https://api.github.com/repos/${owner}/${repo}`);
      if (!response.ok) throw new Error(`GitHub API returned ${response.status}`);
      const json = await response.json();
      const meta = {
        htmlUrl: json.html_url,
        description: json.description,
        language: json.language,
        stars: json.stargazers_count,
        forks: json.forks_count,
        updatedAt: json.updated_at,
        visibility: json.visibility
      };
      state.repoMeta.set(repo, meta);
      return meta;
    } catch (error) {
      state.repoMeta.set(repo, null);
      return null;
    }
  }

  function renderSummary() {
    if (!els.summaryGrid) return;
    const categories = state.data.categories || [];
    const projects = state.data.projects || [];
    const categoryCounts = categories.map((category) => ({
      ...category,
      count: projects.filter((project) => project.category === category.id).length
    }));

    els.summaryGrid.innerHTML = categoryCounts.map((category) => `
      <article class="github-card">
        <div class="icon"><i class="bi bi-diagram-3"></i></div>
        <h3>${escapeHtml(category.label)}</h3>
        <p>${escapeHtml(category.summary)}</p>
        <strong>${category.count} project${category.count === 1 ? '' : 's'}</strong>
      </article>
    `).join('');
  }

  function renderFilters() {
    if (!els.filterRow) return;
    const categories = state.data.categories || [];
    const filters = [{ id: 'all', label: 'All' }, ...categories.map(({ id, label }) => ({ id, label }))];

    els.filterRow.innerHTML = filters.map((filter) => `
      <button class="filter-btn ${filter.id === state.activeFilter ? 'active' : ''}" type="button" data-filter="${escapeHtml(filter.id)}">
        ${escapeHtml(filter.label)}
      </button>
    `).join('');

    els.filterRow.querySelectorAll('[data-filter]').forEach((button) => {
      button.addEventListener('click', () => {
        state.activeFilter = button.dataset.filter;
        renderFilters();
        renderProjects();
      });
    });
  }

  function metaMarkup(project, meta, category) {
    const parts = [
      `<span><i class="bi bi-folder2-open"></i> ${escapeHtml(category?.label || project.category)}</span>`
    ];

    if (meta?.language) {
      parts.push(`<span><i class="bi bi-code-slash"></i> ${escapeHtml(meta.language)}</span>`);
    }

    if (typeof meta?.stars === 'number') {
      parts.push(`<span class="repo-stat"><i class="bi bi-star"></i> ${meta.stars}</span>`);
    }

    if (typeof meta?.forks === 'number') {
      parts.push(`<span class="repo-stat"><i class="bi bi-diagram-2"></i> ${meta.forks}</span>`);
    }

    const updated = formatDate(meta?.updatedAt);
    if (updated) {
      parts.push(`<span><i class="bi bi-clock-history"></i> Updated ${escapeHtml(updated)}</span>`);
    }

    return parts.join('');
  }

  function renderTagList(items = [], className = 'stack-list') {
    return `<ul class="${className}">${items.map((item) => `<li>${escapeHtml(item)}</li>`).join('')}</ul>`;
  }

  function renderProjectCard(project, categories) {
    const meta = state.repoMeta.get(project.repository);
    const category = categories.get(project.category);
    const repoUrl = project.github_url || meta?.htmlUrl || '';
    const liveUrl = project.live_url || '';
    const description = meta?.description || project.tagline;

    return `
      <article class="project-card" data-category="${escapeHtml(project.category)}">
        <div class="project-card-body">
          <div class="project-card-header">
            <h3>${escapeHtml(project.title)}</h3>
            <span class="status-badge status-${escapeHtml(project.status)}">${escapeHtml(statusLabel(project.status))}</span>
          </div>
          <div class="project-meta">${metaMarkup(project, meta, category)}</div>
          <p>${escapeHtml(description)}</p>
          <p><strong>Career signal:</strong> ${escapeHtml(project.why_it_matters)}</p>
          ${renderTagList(project.technical_stack || [])}
          ${project.deliverables?.length ? `<h4 class="mt-3 mb-0 fs-6 text-primary">Evidence to add</h4>${renderTagList(project.deliverables.slice(0, 4), 'deliverable-list')}` : ''}
          <div class="project-actions">
            ${repoUrl ? `<a class="btn btn-outline-primary" href="${escapeHtml(repoUrl)}" target="_blank" rel="noopener"><i class="bi bi-github"></i> Repository</a>` : '<span class="btn btn-outline-secondary disabled">Repository planned</span>'}
            ${liveUrl ? `<a class="btn btn-primary" href="${escapeHtml(liveUrl)}" target="_blank" rel="noopener"><i class="bi bi-box-arrow-up-right"></i> Live demo</a>` : ''}
          </div>
        </div>
      </article>
    `;
  }

  function renderProjects() {
    if (!els.projectGrid) return;
    const categories = categoryById();
    const projects = (state.data.projects || [])
      .filter((project) => state.activeFilter === 'all' || project.category === state.activeFilter)
      .sort((a, b) => a.priority - b.priority);

    els.projectGrid.innerHTML = projects.map((project) => renderProjectCard(project, categories)).join('');
  }

  async function hydrateRepoMetadata() {
    const owner = state.data.profile?.github_owner || OWNER_FALLBACK;
    const projects = state.data.projects || [];

    for (const project of projects) {
      if (project.repository && project.github_url) {
        await getRepoMeta(owner, project.repository);
      }
    }
  }

  async function init() {
    if (!els.projectGrid) return;

    try {
      state.data = await getJson(DATA_URL);
      renderSummary();
      renderFilters();
      renderProjects();
      await hydrateRepoMetadata();
      renderProjects();
      if (els.updated) {
        els.updated.textContent = `Project data loaded from ${DATA_URL}`;
      }
    } catch (error) {
      els.projectGrid.innerHTML = `
        <div class="github-note">
          <h3>Project data could not be loaded</h3>
          <p>Please check <code>${escapeHtml(DATA_URL)}</code> and the browser console.</p>
        </div>
      `;
      console.error(error);
    }
  }

  document.addEventListener('DOMContentLoaded', init);
})();
