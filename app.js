document.addEventListener('DOMContentLoaded', () => {
    // Navigation logic
    const navItems = document.querySelectorAll('.nav-item');
    const views = document.querySelectorAll('.view');

    navItems.forEach(item => {
        item.addEventListener('click', () => {
            const viewId = item.getAttribute('data-view');
            if (!viewId) return;

            // Update active nav
            navItems.forEach(i => i.classList.remove('active'));
            item.classList.add('active');

            // Update visible view
            views.forEach(v => {
                v.classList.remove('active');
                if (v.id === `view-${viewId}`) {
                    v.classList.add('active');
                }
            });
        });
    });

    // Mock Data for Step 1 (M1 & M2)
    const variantsData = [
        { rank: 1, gene: 'PIK3CA', variant: 'H1047R', frequency: 320, relevance: 'High Pathogenic' },
        { rank: 2, gene: 'TP53', variant: 'R175H', frequency: 244, relevance: 'High Pathogenic' },
        { rank: 3, gene: 'PIK3CA', variant: 'E545K', frequency: 201, relevance: 'High Pathogenic' },
        { rank: 4, gene: 'CDH1', variant: 'E243K', frequency: 143, relevance: 'VUS' },
        { rank: 5, gene: 'GATA3', variant: 'P409fs', frequency: 119, relevance: 'Pathogenic' },
        { rank: 6, gene: 'MAP3K1', variant: 'E332K', frequency: 83, relevance: 'VUS' },
        { rank: 7, gene: 'PIK3CA', variant: 'H1047L', frequency: 72, relevance: 'Pathogenic' },
        { rank: 8, gene: 'AKT1', variant: 'E17K', frequency: 54, relevance: 'Pathogenic' },
        { rank: 9, gene: 'PTEN', variant: 'R130Q', frequency: 48, relevance: 'Pathogenic' },
        { rank: 10, gene: 'ESR1', variant: 'D538G', frequency: 39, relevance: 'Pathogenic' }
    ];

    const tableBody = document.getElementById('variants-body');
    if (tableBody) {
        variantsData.forEach(item => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td style="color: var(--text-tertiary)">${item.rank}</td>
                <td><span class="gene-badge">${item.gene}</span></td>
                <td style="font-family: monospace; font-weight: 500;">${item.variant}</td>
                <td>
                    <div style="display: flex; align-items: center; gap: 10px;">
                        <div class="frequency-bar-bg">
                            <div class="frequency-bar-fill" style="width: ${(item.frequency / 320) * 100}%"></div>
                        </div>
                        <span style="font-size: 0.8rem; font-weight: 600;">${item.frequency}</span>
                    </div>
                </td>
                <td>
                    <span style="font-size: 0.75rem; padding: 4px 8px; border-radius: 4px; background: ${item.relevance === 'VUS' ? 'rgba(245, 158, 11, 0.1)' : 'rgba(239, 68, 68, 0.1)'}; color: ${item.relevance === 'VUS' ? '#fbbf24' : '#f87171'}; font-weight: 600;">
                        ${item.relevance}
                    </span>
                </td>
            `;
            tableBody.appendChild(row);
        });
    }

    // Sync button interaction
    const btnSync = document.getElementById('btn-sync');
    if (btnSync) {
        btnSync.addEventListener('click', () => {
            btnSync.innerHTML = '<span>Sincronizando...</span>';
            setTimeout(() => {
                btnSync.innerHTML = '<span>Sincronizado ✓</span>';
                btnSync.style.background = 'var(--success)';
            }, 1500);
        });
    }
});
