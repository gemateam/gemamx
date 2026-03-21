/**
 * analysis.js - Plotly Interactive Visualizations for Mutation Distribution
 */

window.loadAnalysisData = async function() {
    const plotGenes = document.getElementById('plot-genes');
    const plotTypes = document.getElementById('plot-types');

    if (!plotGenes || !plotTypes) return;

    // Show loading state
    plotGenes.innerHTML = '<div style="display:flex; height:100%; align-items:center; justify-content:center; color:var(--text-tertiary)">Cargando analítica...</div>';
    plotTypes.innerHTML = '<div style="display:flex; height:100%; align-items:center; justify-content:center; color:var(--text-tertiary)">Cargando analítica...</div>';

    try {
        const response = await fetch('data/mutation_stats.json');
        if (!response.ok) throw new Error('No se pudo cargar mutation_stats.json');
        
        const stats = await response.json();

        // 1. Gene Distribution Chart (Horizontal Bar)
        const geneLabels = Object.keys(stats.gene_distribution).reverse();
        const geneValues = Object.values(stats.gene_distribution).reverse();

        const geneData = [{
            type: 'bar',
            x: geneValues,
            y: geneLabels,
            orientation: 'h',
            marker: {
                color: 'rgba(13, 115, 119, 0.6)',
                line: {
                    color: 'rgba(13, 115, 119, 1.0)',
                    width: 1
                }
            },
            hovertemplate: '<b>Gen %{y}</b><br>Mutaciones: %{x}<extra></extra>'
        }];

        const geneLayout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#e5e7eb', family: 'Inter, sans-serif' },
            margin: { l: 80, r: 20, t: 10, b: 40 },
            xaxis: { gridcolor: 'rgba(255,255,255,0.05)', title: 'Frecuencia' },
            yaxis: { gridcolor: 'rgba(255,255,255,0.05)' }
        };

        Plotly.newPlot('plot-genes', geneData, geneLayout, { responsive: true, displayModeBar: false });

        // 2. Mutation Type Distribution (Donut Chart)
        const typeLabels = Object.keys(stats.type_distribution);
        const typeValues = Object.values(stats.type_distribution);

        const typeData = [{
            values: typeValues,
            labels: typeLabels,
            type: 'pie',
            hole: 0.4,
            marker: {
                colors: [
                    '#0d7377', '#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6'
                ]
            },
            textinfo: 'percent',
            hoverinfo: 'label+value',
            insidetextorientation: 'radial'
        }];

        const typeLayout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#e5e7eb', family: 'Inter, sans-serif' },
            margin: { l: 20, r: 20, t: 20, b: 20 },
            showlegend: true,
            legend: {
                orientation: 'h',
                x: 0,
                y: -0.1
            }
        };

        if (!typeLabels.length || !typeValues.length) {
            plotTypes.innerHTML = '<p style="color:var(--text-tertiary); text-align:center; padding:40px;">Sin datos de tipo de variante disponibles.</p>';
        } else {
            Plotly.newPlot('plot-types', typeData, typeLayout, { responsive: true, displayModeBar: false });
        }

    } catch (err) {
        console.error("Error al cargar analítica:", err);
        plotGenes.innerHTML = `<p style="color:#ef4444; text-align:center; padding:20px;">Error: ${err.message}</p>`;
        plotTypes.innerHTML = `<p style="color:#ef4444; text-align:center; padding:20px;">Error: No se encontró data/mutation_stats.json. Corre ingest_data.py primero.</p>`;
    }
};

// Also load on first call if view is already active (rare but possible)
if (document.getElementById('view-analysis')?.classList.contains('active')) {
    window.loadAnalysisData();
}
