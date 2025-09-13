// MicrobiotaMR - Interactive Data Explorer
class MicrobiotaMRApp {
    constructor() {
        this.data = [];
        this.filteredData = [];
        this.currentPage = 1;
        this.recordsPerPage = 50;
        this.sortField = null;
        this.sortDirection = 'asc';
        this.chart = null;
        
        this.init();
    }

    async init() {
        await this.loadData();
        this.setupEventListeners();
        this.populateDropdowns();
        this.updateStats();
        this.displayResults();
    }

    async loadData() {
        try {
            const response = await fetch('data.txt');
            const text = await response.text();
            const lines = text.trim().split('\n');
            
            // Parse header
            const headers = lines[0].split('\t');
            console.log('Headers:', headers);
            
            // Parse data
            this.data = lines.slice(1).map((line, index) => {
                const values = line.split('\t');
                return {
                    gene_id: values[0] || '',
                    microbiota: values[1] || '',
                    se: parseFloat(values[3]) || 0,
                    pval: parseFloat(values[2]) || 1,
                    chr: values[6] || '',
                    gene_name: values[4] || '',
                    tissue: values[5] || '',
                    diseases: values[7] || 'NA',
                    id: index
                };
            });
            
            this.filteredData = [...this.data];
            console.log(`Loaded ${this.data.length} records`);
        } catch (error) {
            console.error('Error loading data:', error);
            alert('Error loading data file. Please ensure data.txt is accessible.');
        }
    }

    setupEventListeners() {
        // Search inputs with real-time filtering
        ['gene-search', 'gene-name-search', 'microbiota-search', 'tissue-search'].forEach(id => {
            const input = document.getElementById(id);
            input.addEventListener('input', () => this.handleSearch());
        });

        // Dropdown changes
        ['gene-dropdown', 'gene-name-dropdown', 'microbiota-dropdown', 'tissue-dropdown', 'disease-filter'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => this.handleSearch());
        });

        // NA checkbox
        document.getElementById('include-na').addEventListener('change', () => this.handleSearch());

        // Buttons
        document.getElementById('apply-filters').addEventListener('click', () => this.applyFilters());
        document.getElementById('reset-filters').addEventListener('click', () => this.resetFilters());
        document.getElementById('export-data').addEventListener('click', () => this.exportData());
        document.getElementById('generate-viz').addEventListener('click', () => this.generateVisualization());

        // Table sorting
        document.querySelectorAll('th[data-sort]').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th.dataset.sort));
        });

        // Pagination
        document.getElementById('prev-page').addEventListener('click', () => this.previousPage());
        document.getElementById('next-page').addEventListener('click', () => this.nextPage());
        document.getElementById('records-per-page').addEventListener('change', (e) => {
            this.recordsPerPage = parseInt(e.target.value);
            this.currentPage = 1;
            this.displayResults();
        });
    }

    populateDropdowns() {
        const geneIds = [...new Set(this.data.map(d => d.gene_id))].sort();
        const geneNames = [...new Set(this.data.map(d => d.gene_name))].sort();
        const microbiotas = [...new Set(this.data.map(d => d.microbiota))].sort();
        const tissues = [...new Set(this.data.map(d => d.tissue))].sort();
        const diseases = [...new Set(this.data.map(d => d.diseases).filter(d => d !== 'NA'))].sort();

        this.populateDropdown('gene-dropdown', geneIds);
        this.populateDropdown('gene-name-dropdown', geneNames);
        this.populateDropdown('microbiota-dropdown', microbiotas);
        this.populateDropdown('tissue-dropdown', tissues);
        this.populateDropdown('disease-filter', diseases);
    }

    populateDropdown(elementId, options) {
        const select = document.getElementById(elementId);
        const defaultOption = select.querySelector('option[value=""]');
        select.innerHTML = '';
        select.appendChild(defaultOption);
        
        options.forEach(option => {
            if (option && option.trim()) {
                const optionElement = document.createElement('option');
                optionElement.value = option;
                optionElement.textContent = option;
                select.appendChild(optionElement);
            }
        });
    }

    handleSearch() {
        const filters = this.getFilterValues();
        this.applyFiltersWithValues(filters);
    }

    getFilterValues() {
        return {
            geneId: document.getElementById('gene-search').value.toLowerCase().trim(),
            geneName: document.getElementById('gene-name-search').value.toLowerCase().trim(),
            microbiota: document.getElementById('microbiota-search').value.toLowerCase().trim(),
            tissue: document.getElementById('tissue-search').value.toLowerCase().trim(),
            geneDropdown: document.getElementById('gene-dropdown').value,
            geneNameDropdown: document.getElementById('gene-name-dropdown').value,
            microbiotaDropdown: document.getElementById('microbiota-dropdown').value,
            tissueDropdown: document.getElementById('tissue-dropdown').value,
            diseaseFilter: document.getElementById('disease-filter').value,
            includeNA: document.getElementById('include-na').checked
        };
    }

    applyFiltersWithValues(filters) {
        this.filteredData = this.data.filter(row => {
            // Text search filters
            const geneIdMatch = !filters.geneId || row.gene_id.toLowerCase().includes(filters.geneId);
            const geneNameMatch = !filters.geneName || row.gene_name.toLowerCase().includes(filters.geneName);
            const microbiotaMatch = !filters.microbiota || row.microbiota.toLowerCase().includes(filters.microbiota);
            const tissueMatch = !filters.tissue || row.tissue.toLowerCase().includes(filters.tissue);

            // Dropdown filters
            const geneDropdownMatch = !filters.geneDropdown || row.gene_id === filters.geneDropdown;
            const geneNameDropdownMatch = !filters.geneNameDropdown || row.gene_name === filters.geneNameDropdown;
            const microbiotaDropdownMatch = !filters.microbiotaDropdown || row.microbiota === filters.microbiotaDropdown;
            const tissueDropdownMatch = !filters.tissueDropdown || row.tissue === filters.tissueDropdown;
            
            // Disease filter logic
            let diseaseMatch = true;
            if (filters.diseaseFilter) {
                diseaseMatch = row.diseases === filters.diseaseFilter;
            } else {
                // If no specific disease selected, check NA inclusion
                if (!filters.includeNA && row.diseases === 'NA') {
                    diseaseMatch = false;
                }
            }

            return geneIdMatch && geneNameMatch && microbiotaMatch && tissueMatch &&
                   geneDropdownMatch && geneNameDropdownMatch && microbiotaDropdownMatch &&
                   tissueDropdownMatch && diseaseMatch;
        });

        // Update dependent dropdowns based on current filters
        this.updateDependentDropdowns(filters);
        
        this.currentPage = 1;
        this.updateStats();
        this.displayResults();
    }

    updateDependentDropdowns(currentFilters) {
        // Get available options based on current filtered data
        const availableGeneIds = [...new Set(this.filteredData.map(d => d.gene_id))].sort();
        const availableGeneNames = [...new Set(this.filteredData.map(d => d.gene_name))].sort();
        const availableMicrobiome = [...new Set(this.filteredData.map(d => d.microbiota))].sort();
        const availableTissues = [...new Set(this.filteredData.map(d => d.tissue))].sort();
        const availableDiseases = [...new Set(this.filteredData.map(d => d.diseases).filter(d => d !== 'NA'))].sort();

        // Update dropdowns while preserving selected values
        this.updateDropdownOptions('gene-dropdown', availableGeneIds, currentFilters.geneDropdown);
        this.updateDropdownOptions('gene-name-dropdown', availableGeneNames, currentFilters.geneNameDropdown);
        this.updateDropdownOptions('microbiota-dropdown', availableMicrobiome, currentFilters.microbiotaDropdown);
        this.updateDropdownOptions('tissue-dropdown', availableTissues, currentFilters.tissueDropdown);
        this.updateDropdownOptions('disease-filter', availableDiseases, currentFilters.diseaseFilter);
    }

    updateDropdownOptions(elementId, options, selectedValue) {
        const select = document.getElementById(elementId);
        const defaultText = select.querySelector('option[value=""]').textContent;
        
        select.innerHTML = `<option value="">${defaultText}</option>`;
        
        options.forEach(option => {
            if (option && option.trim()) {
                const optionElement = document.createElement('option');
                optionElement.value = option;
                optionElement.textContent = option;
                if (option === selectedValue) {
                    optionElement.selected = true;
                }
                select.appendChild(optionElement);
            }
        });
    }

    applyFilters() {
        const filters = this.getFilterValues();
        this.applyFiltersWithValues(filters);
    }

    resetFilters() {
        // Clear all inputs
        ['gene-search', 'gene-name-search', 'microbiota-search', 'tissue-search'].forEach(id => {
            document.getElementById(id).value = '';
        });

        // Reset all dropdowns
        ['gene-dropdown', 'gene-name-dropdown', 'microbiota-dropdown', 'tissue-dropdown', 'disease-filter'].forEach(id => {
            document.getElementById(id).value = '';
        });

        // Reset NA checkbox to checked
        document.getElementById('include-na').checked = true;

        // Reset filtered data
        this.filteredData = [...this.data];
        this.currentPage = 1;
        
        // Repopulate all dropdowns
        this.populateDropdowns();
        this.updateStats();
        this.displayResults();
    }

    sortTable(field) {
        if (this.sortField === field) {
            this.sortDirection = this.sortDirection === 'asc' ? 'desc' : 'asc';
        } else {
            this.sortField = field;
            this.sortDirection = 'asc';
        }

        this.filteredData.sort((a, b) => {
            let aVal = a[field];
            let bVal = b[field];

            // Handle numeric fields
            if (field === 'pval' || field === 'se' || field === 'chr') {
                aVal = parseFloat(aVal) || 0;
                bVal = parseFloat(bVal) || 0;
            } else {
                aVal = String(aVal).toLowerCase();
                bVal = String(bVal).toLowerCase();
            }

            if (aVal < bVal) return this.sortDirection === 'asc' ? -1 : 1;
            if (aVal > bVal) return this.sortDirection === 'asc' ? 1 : -1;
            return 0;
        });

        this.updateSortIcons();
        this.displayResults();
    }

    updateSortIcons() {
        document.querySelectorAll('.sort-icon').forEach(icon => {
            icon.textContent = '⇅';
        });

        if (this.sortField) {
            const th = document.querySelector(`th[data-sort="${this.sortField}"] .sort-icon`);
            if (th) {
                th.textContent = this.sortDirection === 'asc' ? '▲' : '▼';
            }
        }
    }

    displayResults() {
        const tbody = document.getElementById('results-tbody');
        tbody.innerHTML = '';

        const startIndex = (this.currentPage - 1) * this.recordsPerPage;
        const endIndex = Math.min(startIndex + this.recordsPerPage, this.filteredData.length);

        for (let i = startIndex; i < endIndex; i++) {
            const row = this.filteredData[i];
            const tr = document.createElement('tr');

            tr.innerHTML = `
                <td>${row.gene_id}</td>
                <td>${row.gene_name}</td>
                <td title="${row.microbiota}">${this.truncateText(row.microbiota, 30)}</td>
                <td>${row.tissue}</td>
                <td>${row.diseases}</td>
                <td>${row.chr}</td>
                <td>${row.pval.toExponential(3)}</td>
                <td>${row.se.toFixed(6)}</td>
            `;

            tbody.appendChild(tr);
        }

        this.updatePagination();
    }

    truncateText(text, maxLength) {
        return text.length > maxLength ? text.substring(0, maxLength) + '...' : text;
    }

    updatePagination() {
        const totalPages = Math.ceil(this.filteredData.length / this.recordsPerPage);
        document.getElementById('page-info').textContent = `Page ${this.currentPage} of ${totalPages}`;
        
        document.getElementById('prev-page').disabled = this.currentPage === 1;
        document.getElementById('next-page').disabled = this.currentPage === totalPages || totalPages === 0;
    }

    previousPage() {
        if (this.currentPage > 1) {
            this.currentPage--;
            this.displayResults();
        }
    }

    nextPage() {
        const totalPages = Math.ceil(this.filteredData.length / this.recordsPerPage);
        if (this.currentPage < totalPages) {
            this.currentPage++;
            this.displayResults();
        }
    }

    updateStats() {
        const uniqueGenes = new Set(this.filteredData.map(d => d.gene_id)).size;
        const uniqueMicrobiome = new Set(this.filteredData.map(d => d.microbiota)).size;

        document.getElementById('total-records').textContent = this.data.length.toLocaleString();
        document.getElementById('filtered-records').textContent = this.filteredData.length.toLocaleString();
        document.getElementById('unique-genes').textContent = uniqueGenes.toLocaleString();
        document.getElementById('unique-microbiota').textContent = uniqueMicrobiome.toLocaleString();
    }

    generateVisualization() {
        const vizType = document.getElementById('viz-type').value;
        
        if (this.filteredData.length === 0) {
            alert('No data to visualize. Please adjust your filters.');
            return;
        }

        switch (vizType) {
            case 'scatter':
                this.createScatterPlot();
                break;
            case 'heatmap':
                this.createHeatmap();
                break;
        }
    }

    createScatterPlot() {
        // Clear D3 container and show canvas
        const d3Container = document.getElementById('d3-container');
        d3Container.innerHTML = '';
        d3Container.style.display = 'none';
        
        const canvas = document.getElementById('main-chart');
        canvas.style.display = 'block';
        const ctx = canvas.getContext('2d');
        
        if (this.chart) {
            this.chart.destroy();
        }

        const data = this.filteredData.map(d => ({
            x: Math.log10(d.se),
            y: -Math.log10(d.pval),
            gene: d.gene_name,
            microbiota: d.microbiota
        }));

        this.chart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Gene-microbiota Associations',
                    data: data,
                    backgroundColor: 'rgba(44, 90, 160, 0.6)',
                    borderColor: 'rgba(44, 90, 160, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'log10(Standard Error)'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: '-log10(P-value)'
                        }
                    }
                },
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const point = context.raw;
                                return [
                                    `Gene: ${point.gene}`,
                                    `Microbiota: ${point.microbiota}`,
                                    `SE: ${Math.pow(10, point.x).toExponential(3)}`,
                                    `P-value: ${Math.pow(10, -point.y).toExponential(3)}`
                                ];
                            }
                        }
                    }
                }
            }
        });
    }

    createHeatmap() {
        // Destroy any existing Chart.js chart
        if (this.chart) {
            this.chart.destroy();
            this.chart = null;
        }
        
        // Clear chart canvas and use D3 for heatmap
        const canvas = document.getElementById('main-chart');
        canvas.style.display = 'none';
        
        const container = document.getElementById('d3-container');
        container.innerHTML = '';
        container.style.display = 'block';

        // Prepare data for heatmap
        const topGenes = [...new Set(this.filteredData.map(d => d.gene_name))].slice(0, 20);
        const topMicrobiota = [...new Set(this.filteredData.map(d => d.microbiota))].slice(0, 15);

        const heatmapData = [];
        topGenes.forEach((gene, geneIdx) => {
            topMicrobiota.forEach((microbiota, microIdx) => {
                const records = this.filteredData.filter(d => 
                    d.gene_name === gene && d.microbiota === microbiota
                );
                if (records.length > 0) {
                    const avgPval = records.reduce((sum, r) => sum + r.pval, 0) / records.length;
                    heatmapData.push({
                        gene: gene,
                        microbiota: microbiota,
                        geneIdx: geneIdx,
                        microIdx: microIdx,
                        pval: avgPval,
                        logPval: -Math.log10(avgPval)
                    });
                }
            });
        });

        this.renderD3Heatmap(container, heatmapData, topGenes, topMicrobiota);
    }

    renderD3Heatmap(container, data, genes, microbiota) {
        const margin = {top: 80, right: 50, bottom: 200, left: 150};
        const width = 800 - margin.left - margin.right;
        const height = 600 - margin.top - margin.bottom;

        const svg = d3.select(container)
            .append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom);

        const g = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        const xScale = d3.scaleBand()
            .range([0, width])
            .domain(microbiota)
            .padding(0.1);

        const yScale = d3.scaleBand()
            .range([height, 0])
            .domain(genes)
            .padding(0.1);

        const colorScale = d3.scaleSequential(d3.interpolateBlues)
            .domain([0, d3.max(data, d => d.logPval)]);

        // Add rectangles
        g.selectAll('rect')
            .data(data)
            .enter()
            .append('rect')
            .attr('x', d => xScale(d.microbiota))
            .attr('y', d => yScale(d.gene))
            .attr('width', xScale.bandwidth())
            .attr('height', yScale.bandwidth())
            .style('fill', d => colorScale(d.logPval))
            .append('title')
            .text(d => `${d.gene} - ${d.microbiota}\nP-value: ${d.pval.toExponential(3)}`);

        // Add x-axis
        g.append('g')
            .attr('transform', `translate(0,${height})`)
            .call(d3.axisBottom(xScale))
            .selectAll('text')
            .style('text-anchor', 'end')
            .attr('dx', '-.8em')
            .attr('dy', '.15em')
            .attr('transform', 'rotate(-45)');

        // Add y-axis
        g.append('g')
            .call(d3.axisLeft(yScale));

        // Add title
        svg.append('text')
            .attr('x', (width + margin.left + margin.right) / 2)
            .attr('y', margin.top / 2)
            .attr('text-anchor', 'middle')
            .style('font-size', '16px')
            .style('font-weight', 'bold')
            .text('Gene-Microbiota Association Heatmap');
    }

    exportData() {
        if (this.filteredData.length === 0) {
            alert('No data to export. Please adjust your filters.');
            return;
        }

        const csvHeader = 'Gene ID,Gene Name,Microbiota,Tissue,Disease,Chromosome,P-value,Standard Error\n';
        const csvContent = this.filteredData.map(row => [
            row.gene_id,
            row.gene_name,
            row.microbiota,
            row.tissue,
            row.diseases,
            row.chr,
            row.pval,
            row.se
        ].join(',')).join('\n');

        const blob = new Blob([csvHeader + csvContent], { type: 'text/csv' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `MicrobiotaMR_filtered_${new Date().toISOString().split('T')[0]}.csv`;
        a.click();
        window.URL.revokeObjectURL(url);
    }
}

// Initialize the application when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new MicrobiotaMRApp();
});
