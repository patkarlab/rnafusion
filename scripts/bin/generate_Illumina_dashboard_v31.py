#!/usr/bin/env python3
"""
Gene Fusion Analysis Dashboard Generator

Usage:
    python generate_illumina_dashboard_v31.py <excel_file> <cytoband_file> <output_html>

Example:
    python generate_illumina_dashboard_v31.py sample.xlsx cytoBand.txt sample_dashboard.html
"""

import sys
import json
import pandas as pd
import argparse
from pathlib import Path

def parse_coverage_sheet(excel_file):
    """Parse coverage data from the first sheet"""
    df = pd.read_excel(excel_file, sheet_name=0, header=None)
    coverage_data = []
    
    for _, row in df.iterrows():
        if len(row) >= 6 and pd.notna(row[0]):
            coverage_data.append({
                'chr': str(row[0]),
                'start': int(row[1]) if pd.notna(row[1]) else 0,
                'end': int(row[2]) if pd.notna(row[2]) else 0,
                'target': str(row[3]) if pd.notna(row[3]) else '',
                'strand': str(row[4]) if pd.notna(row[4]) else '',
                'coverage': int(row[5]) if pd.notna(row[5]) else 0
            })
    
    return coverage_data


def extract_gene_name(row, keys):
    """Extract gene name from row using list of possible column keys"""
    for key in keys:
        if key in row and pd.notna(row[key]):
            value = str(row[key]).strip()
            if value not in ['.', '']:
                if '--' in value:
                    return value.split('--')[0]
                return value
    return ''


def extract_gene2_name(row, keys, gene1_keys):
    """Extract gene2 name from row"""
    # Try direct gene2 columns first
    for key in keys:
        if key in row and pd.notna(row[key]):
            value = str(row[key]).strip()
            if value not in ['.', '']:
                return value
    
    # Try fusion format from gene1 columns
    for key in gene1_keys:
        if key in row and pd.notna(row[key]):
            value = str(row[key]).strip()
            if value not in ['.', ''] and '--' in value:
                return value.split('--')[1]
    
    return ''


def extract_value(row, keys, default=''):
    """Extract value from row using list of possible keys"""
    for key in keys:
        if key in row and pd.notna(row[key]):
            return str(row[key])
    return default


def extract_number(row, keys, default=0):
    """Extract numeric value from row"""
    for key in keys:
        if key in row and pd.notna(row[key]):
            try:
                return int(float(row[key]))
            except (ValueError, TypeError):
                pass
    return default


def extract_chr_num(breakpoint_str):
    """Extract chromosome number from breakpoint string for comparison"""
    try:
        chr_part = breakpoint_str.split(':')[0]
        chr_part = chr_part.replace('chr', '')
        # Handle X, Y chromosomes
        if chr_part == 'X':
            return 23
        elif chr_part == 'Y':
            return 24
        return int(chr_part)
    except:
        return 999  # Put invalid ones at end


def parse_fusion_callers(fusion_dfs):
    """Parse fusion data from already-loaded dataframe dict"""

    fusion_map = {}

    for caller, df in fusion_dfs.items():
        print(f"Processing {caller} sheet with {len(df)} rows")

        for _, row in df.iterrows():

            # ===================================
            # ARRIBA
            # ===================================
            if caller == "Arriba":
                gene1 = str(row.get("#gene1", "")).strip()
                gene2 = str(row.get("gene2", "")).strip()

                breakpoint1_str = str(row.get("breakpoint1", "")).strip()
                breakpoint2_str = str(row.get("breakpoint2", "")).strip()

                split_reads = extract_number(row, ["split_reads1"])
                discordant = extract_number(row, ["discordant_mates"])

                confidence = row.get("confidence", "")
                retained_domains = row.get("retained_protein_domains", "")
                coding_frames = row.get("reading_frame", "")
                annotation = ""

            # ===================================
            # SQUID
            # ===================================
            elif caller == "Squid":
                # Gene names
                fused = str(row.get("FusedGenes", "")).strip()
                if ":" not in fused:
                    continue

                gene1, gene2 = [x.strip() for x in fused.split(":", 1)]

                # Breakpoints
                chrom1 = str(row.get("# chrom1", "")).replace("chr", "").strip()
                chrom2 = str(row.get("chrom2", "")).replace("chr", "").strip()
                start1 = str(row.get("start1", "")).strip()
                start2 = str(row.get("start2", "")).strip()

                if not (chrom1 and chrom2 and start1 and start2):
                    continue

                bp1_raw = f"chr{chrom1}:{start1}"
                bp2_raw = f"chr{chrom2}:{start2}"

                # Order
                if extract_chr_num(bp1_raw) <= extract_chr_num(bp2_raw):
                    breakpoint1_str = bp1_raw
                    breakpoint2_str = bp2_raw
                else:
                    breakpoint1_str = bp2_raw
                    breakpoint2_str = bp1_raw
                    gene1, gene2 = gene2, gene1

                split_reads = 0
                discordant = extract_number(row, ["score"])
                confidence = ""
                retained_domains = ""
                coding_frames = ""
                annotation = ""

            # ===================================
            # PIZZLY
            # ===================================
            elif caller == "Pizzly":
                gene1 = str(row.get("geneA.name", "")).strip()
                gene2 = str(row.get("geneB.name", "")).strip()

                split_reads = extract_number(row, ["splitcount"])
                discordant = extract_number(row, ["paircount"])

                breakpoint1_str = ""
                breakpoint2_str = ""

                confidence = ""
                retained_domains = ""
                coding_frames = ""
                annotation = ""

            # ===================================
            # FUSIONCATCHER
            # ===================================
            elif caller == "FusionCatcher":
                gene1 = str(row.get("Gene_1_symbol(5end_fusion_partner)", "")).strip()
                gene2 = str(row.get("Gene_2_symbol(3end_fusion_partner)", "")).strip()

                breakpoint1_str = str(row.get("Fusion_point_for_gene_1(5end_fusion_partner)", "")).strip()
                breakpoint2_str = str(row.get("Fusion_point_for_gene_2(3end_fusion_partner)", "")).strip()

                split_reads = extract_number(row, ["Spanning_pairs"])
                discordant = extract_number(row, ["Counts_of_common_mapping_reads"])

                annotation = row.get("Fusion_description", "")
                confidence = ""
                retained_domains = ""
                coding_frames = ""

            # ===================================
            # STAR-FUSION
            # ===================================
            elif caller == "STAR-Fusion":
                fusion_field = str(row.get("#FusionName", "")).strip()
                if "--" in fusion_field:
                    g1, g2 = fusion_field.split("--", 1)
                    gene1 = g1.strip()
                    gene2 = g2.strip()
                else:
                    continue

                breakpoint1_str = str(row.get("LeftBreakpoint", "")).strip()
                breakpoint2_str = str(row.get("RightBreakpoint", "")).strip()

                split_reads = extract_number(row, ["JunctionReadCount"])
                discordant = extract_number(row, ["SpanningFragCount"])

                annotation = row.get("annots", "")
                confidence = ""
                retained_domains = ""
                coding_frames = ""

            else:
                continue  # unknown caller

            # ===================================
            # Skip if missing gene names
            # ===================================
            if not gene1 or not gene2:
                continue

            # ===================================
            # Normalize key
            # ===================================
            normalized_key = "--".join(sorted([gene1.upper(), gene2.upper()]))
            fusion_display = f"{gene1}--{gene2}"

            if normalized_key not in fusion_map:
                fusion_map[normalized_key] = {
                    "fusionName": fusion_display,
                    "gene1": gene1,
                    "gene2": gene2,
                    "callers": [],
                    "totalReads": 0,
                    "breakpoints": []
                }

            fusion = fusion_map[normalized_key]

            # Add breakpoint entry
            bp = {
                "caller": caller,
                "gene1": gene1,
                "gene2": gene2,
                "breakpoint1": breakpoint1_str,
                "breakpoint2": breakpoint2_str,
                "splitReads": split_reads,
                "discordantMates": discordant,
                "confidence": confidence,
                "annotation": annotation,
                "codingFrames": coding_frames,
                "retainedProteinDomains": retained_domains,
                "type": row.get("type", "")
            }

            # ---------------------------------------
            #   DEDUP: keep only highest-read entry
            # ---------------------------------------
            total_reads = split_reads + discordant
            
            # Check if this caller already reported this fusion
            existing_idx = None
            for idx, existing in enumerate(fusion["breakpoints"]):
                if existing["caller"] == caller:
                    existing_idx = idx
                    break
            
            if existing_idx is not None:
                # Compare read counts
                old_bp = fusion["breakpoints"][existing_idx]
                old_reads = (old_bp["splitReads"] or 0) + (old_bp["discordantMates"] or 0)
            
                if total_reads > old_reads:
                    # Replace with better (higher support)
                    fusion["breakpoints"][existing_idx] = bp
            
                    # Adjust totalReads
                    fusion["totalReads"] -= old_reads
                    fusion["totalReads"] += total_reads
            else:
                # First time this caller reports this fusion
                fusion["breakpoints"].append(bp)
                fusion["callers"].append(caller)
                fusion["totalReads"] += total_reads


    print(f"Total unique fusions: {len(fusion_map)}")
    return list(fusion_map.values())


def load_cytoband_data(cytoband_file):
    valid_chrs = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    cytobands = []
    with open(cytoband_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            chrom = parts[0]

            # KEEP ONLY standard chromosomes
            if chrom not in valid_chrs:
                continue

            cytobands.append({
                'chr': chrom,
                'start': int(parts[1]),
                'end': int(parts[2]),
                'band': parts[3],
                'stain': parts[4]
            })

    return cytobands



def generate_html(coverage_data, fusion_data, cytoband_data, output_file, sample_name):
    """Generate standalone HTML file with embedded data and visualization"""
    
    html_template = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Fusion Analysis Dashboard</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        :root {
            --bg-primary: hsl(200, 80%, 70%);
            --bg-card: hsl(0, 0%, 100%);
            --text-primary: hsl(200, 30%, 10%);
            --text-muted: hsl(200, 15%, 45%);
            --border: hsl(214, 32%, 91%);
            --primary: hsl(217, 91%, 60%);
            --primary-fg: hsl(0, 0%, 100%);
            --accent: hsl(142, 76%, 45%);
            --destructive: hsl(0, 84%, 60%);
            --muted-bg: hsl(210, 40%, 96%);
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: var(--text-primary);
            font-size: 18px;
            font-weight: 600;
            line-height: 1.6;
        }
        
        .container {
	    width: 100%;
	    max-width: 1600px;
	    margin: 0 auto;
	    padding: 2rem 0.5rem;
	}
        
        .header-container {
            background: white;
            border-radius: 0.75rem;
            padding: 1.5rem 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
        }

        
        h1 {
            font-size: 2.75rem;
            font-weight: bold;
            margin-bottom: 0.5rem;
        }
        
        .subtitle {
            color: var(--text-muted);
            font-size: 1.125rem;
            font-weight: 600;
            margin-bottom: 2rem;
        }
        
        .card {
            background: var(--bg-card);
            border-radius: 0.5rem;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        .card-title {
            font-size: 1.375rem;
            font-weight: 700;
            margin-bottom: 1rem;
        }
        
        .filters {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
            margin-bottom: 1rem;
        }
        
        .filter-group label {
            display: block;
            font-size: 1.0625rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
        }
        
        select, input {
	    width: 100%;
	    padding: 0.5rem;

	    border: 1px solid hsl(217, 91%, 40%);   /* dark blue border */
	    border-radius: 0.375rem;

	    background: hsl(210, 40%, 98%);         /* soft raised look background */
	    font-size: 1.0625rem;
	    font-weight: 600;

	    /* NEW: 3D EFFECT */
	    box-shadow: 0 2px 4px rgba(0,0,0,0.2);
	    transition: all 0.2s ease;
	}

	select:hover, input:hover {
	    box-shadow: 0 3px 6px rgba(0,0,0,0.25);
	    background: white;
	}

	select:focus, input:focus {
	    outline: none;
	    box-shadow: 0 4px 10px rgba(0,0,0,0.3);
	    background: white;
	    border-color: var(--primary);
	}

        
        .tabs {
            margin-bottom: 1rem;
        }
        
        .tab-list {
            display: flex;
            gap: 0.5rem;
            border-bottom: 1px solid var(--border);
            margin-bottom: 1rem;
        }
        
        .tab-button {
	    padding: 0.75rem 1.5rem;
	    border: 1px solid rgba(255,255,255,0.2);
	    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
	    cursor: pointer;
	    font-size: 1rem;
	    font-weight: 700;
	    color: white;
	    border-bottom: none;
	    border-radius: 0.5rem 0.5rem 0 0;

	    /* NEW: 3D effect */
	    box-shadow: 0 3px 6px rgba(0,0,0,0.25);
	    transition: all 0.2s ease;
	}

	.tab-button:hover {
	    background: linear-gradient(135deg, #5a67d8 0%, #6b3d94 100%);
	    transform: translateY(-2px);     /* lift on hover */
	    box-shadow: 0 5px 10px rgba(0,0,0,0.3);
	}

	.tab-button.active {
	    background: linear-gradient(135deg, #7c8ff0 0%, #8b5fc7 100%);
	    color: white;

	    border-bottom: 2px solid rgba(255,255,255,0.3);
	    border-top: 2px solid rgba(255,255,255,0.3);
	    font-weight: 600;

	    /* Brighter, more raised */
	    box-shadow: 0 4px 12px rgba(0,0,0,0.35);
	    transform: translateY(-2px);
	}
	
	/* FIRST HEADER ROW — Fusion | Callers | Arriba | ... | Total Reads */
	table thead tr:first-child th {
	    position: sticky;
	    top: 0;
	    z-index: 10;
	    background: #667eea;
	    color: white;
	    padding: 0.75rem;
	    font-weight: 700;
	    text-align: left;

	    box-shadow: 0 3px 6px rgba(0,0,0,0.25);
	    border-bottom: 2px solid #5a67d8;

	    border-top-left-radius: 0.5rem;
	    border-top-right-radius: 0.5rem;

	    transition: all 0.2s ease;
	}

	/* SECOND HEADER ROW — Split | Spanning | ... */
	table thead tr:nth-child(2) th {
	    position: sticky;
	    top: 3.2rem;                  /* adjust if needed */
	    z-index: 9;
	    background: hsl(210, 40%, 92%);
	    color: #333;
	    padding: 0.5rem;
	    font-weight: 700;
	    text-align: left;

	    box-shadow: none;
	    border: 1px solid var(--border);
	}
 

        .tab-content {
            display: none;
        }
        
        .tab-content.active {
            display: block;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            font-size: 1.0625rem;
            font-weight: 600;
        }
        
        th {
            background: #667eea;
            color: white;
            padding: 0.75rem;
            text-align: left;
            font-weight: 700;
            border: 1px solid rgba(255,255,255,0.2);
        }
        
        td {
            padding: 0.75rem;
            border: 1px solid var(--border);
        }
        
        tr:hover {
            background: var(--muted-bg);
        }
        
        .table-container {
	    max-height: 600px;                 /* allows vertical scrolling */
	    overflow-y: auto;                  /* enable vertical scrolling */
	    overflow-x: auto;                  /* keep your horizontal scroll */
	    border-radius: 0.5rem;
	    border: 1px solid var(--border);
	    position: relative; /* required for sticky */
	}

        
        .badge {
            display: inline-block;
            padding: 0.25rem 0.75rem;
            border-radius: 9999px;
            font-size: 0.9375rem;
            font-weight: 700;
            background: var(--primary);
            color: var(--primary-fg);
        }
        
        .badge.secondary {
            background: var(--muted-bg);
            color: var(--text-primary);
        }
        
        .fusion-row {
            cursor: pointer;
        }
        
        .fusion-row:hover {
            background: var(--muted-bg) !important;
        }
        
        .fusion-detail {
            display: none;
            margin-top: 1.5rem;
        }
        
        .fusion-detail.active {
            display: block;
        }
        
        .circos-layout {
            display: flex;
            gap: 1.5rem;
            align-items: flex-start;
        }
        
        .circos-sidebar {
            display: flex;
            flex-direction: column;
            gap: 1rem;
            width: 12rem;
            flex-shrink: 0;
        }
        
        .circos-main {
            flex: 1;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 1rem;
        }
        
        .stat-card {
            background: var(--muted-bg);
            padding: 1rem;
            border-radius: 0.5rem;
            text-align: center;
        }
        
        .stat-value {
            font-size: 2rem;
            font-weight: bold;
            color: var(--primary);
        }
        
        .summary-cards {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1.5rem;
            margin-bottom: 1.5rem;
        }
        
        .summary-card {
            background: var(--bg-card);
            padding: 1.5rem;
            border-radius: 0.75rem;
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1), 0 2px 4px -1px rgba(0,0,0,0.06);
            text-align: center;
            border: 1px solid var(--border);
        }
        
        .summary-card-value {
            font-size: 3rem;
            font-weight: bold;
            color: var(--primary);
            margin-bottom: 0.5rem;
        }
        
        .summary-card-label {
            font-size: 1.125rem;
            font-weight: 700;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }
        
        .stat-label {
            font-size: 1.0625rem;
            font-weight: 700;
            color: var(--text-muted);
        }
        
        #circosPlot {
            display: flex;
            justify-content: center;
            align-items: center;
            background: var(--bg-card);
            border-radius: 0.5rem;
            padding: 1rem;
        }
        
        .text-center { text-align: center; }
        .text-right { text-align: right; }
        .font-bold { font-weight: bold; }
        .text-primary { color: var(--primary); }
        
        .loading {
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 200px;
        }
        
        .spinner {
            border: 3px solid var(--border);
            border-top-color: var(--primary);
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }
        
        @keyframes spin {
            to { transform: rotate(360deg); }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header-container">
            <h1>Gene Fusion Analysis Dashboard</h1>
            <p class="subtitle">Sample: ''' + sample_name + '''</p>
        </div>
        
        <div class="card">
            <div class="card-title">Filters & Sorting</div>
            <div class="filters">
                <div class="filter-group">
                    <label>Minimum Callers</label>
                    <select id="minCallers">
                        <option value="1">1+ callers</option>
                        <option value="2">2+ callers</option>
                        <option value="3">3+ callers</option>
                        <option value="4">4+ callers</option>
                        <option value="5">5+ callers</option>
                    </select>
                </div>
                <div class="filter-group">
                    <label>Minimum Reads</label>
                    <input type="number" id="minReads" value="0" min="0" step="1">
                </div>
                <div class="filter-group">
                    <label>Sort By</label>
                    <select id="sortBy">
                        <option value="callers">Number of Callers</option>
                        <option value="reads">Total Reads</option>
                        <option value="name">Fusion Name</option>
                    </select>
                </div>
            </div>
        </div>
        
        <div class="summary-cards" id="summaryCards">
            <!-- Summary cards will be populated by JavaScript -->
        </div>
        
        <div class="tabs">
            <div class="tab-list">
                <button class="tab-button active" data-tab="fusions">Fusion Results</button>
                <button class="tab-button" data-tab="coverage">Coverage Data</button>
            </div>
            
            <div id="fusions-tab" class="tab-content active">
                <div class="card">
                    <div id="fusionsTable"></div>
                </div>
            
                <div id="fusionDetail" class="fusion-detail">
            
                    <!-- FUSION DETAILS (Title + Stats + Callers) -->
                    <div class="card">
                        <div class="card-title" id="detailTitle"></div>
            
                        <!-- Stats -->
                        <div id="detailStats"></div>
            
                        <!-- Callers list -->
                        <p class="text-center text-muted" 
                           id="callersList" 
                           style="font-size: 0.875rem; margin-top: 0.75rem;">
                        </p>
            
                        <!-- BREAKPOINT TABLE (MOVED HERE, BEFORE CIRCOS) -->
                        <div class="card" style="margin-top: 1rem;">
                            <div class="card-title">Breakpoints</div>
                            <div id="breakpointsTable"></div>
                        </div>
            
                        <!-- FULL CIRCOS PANEL (now rendered AFTER breakpoints) -->
                        <div class="circos-layout" style="margin-top: 1.5rem;">
                            <div class="circos-sidebar">
                                <!-- left empty intentionally or can contain stats later -->
                            </div>
            
                            <div class="circos-main">
                                <div id="circosPlot"></div>
            
                                <!-- keep fusion gene labels & "detected by" inside circos block if needed -->
                            </div>
                        </div>
            
                    </div> <!-- end fusion detail card -->
            
                </div>
            </div>

            
            <div id="coverage-tab" class="tab-content">
                <div class="card">
                    <div id="coverageTable"></div>
                </div>
            </div>
        </div>
    </div>
    
    
    <script>
        // Embedded data
        const coverageData = ''' + json.dumps(coverage_data) + ''';
        const fusionData = ''' + json.dumps(fusion_data) + ''';
        const cytobandData = ''' + json.dumps(cytoband_data) + ''';
        
        let currentFusion = null;
        
        // Tab switching
        document.querySelectorAll('.tab-button').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const tab = e.target.dataset.tab;
                const parent = e.target.closest('.tabs');
                
                parent.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
                parent.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                
                e.target.classList.add('active');
                parent.querySelector(`#${tab}-tab`).classList.add('active');
            });
        });
        
        // Calculate summary statistics
        function updateSummaryCards() {
            const totalFusions = fusionData.length;
            const highConfidenceFusions = fusionData.filter(f => f.totalReads > 10).length;
            
            document.getElementById('summaryCards').innerHTML = `
                <div class="summary-card">
                    <div class="summary-card-value">${totalFusions}</div>
                    <div class="summary-card-label">Total Fusions</div>
                </div>
                <div class="summary-card">
                    <div class="summary-card-value">${highConfidenceFusions}</div>
                    <div class="summary-card-label">High Confidence (>10 reads)</div>
                </div>
            `;
        }
        
        // Filter and sort
        function getFilteredFusions() {
            const minCallers = parseInt(document.getElementById('minCallers').value);
            const minReads = parseInt(document.getElementById('minReads').value);
            const sortBy = document.getElementById('sortBy').value;
            
            let filtered = fusionData.filter(f => 
                f.callers.length >= minCallers && f.totalReads >= minReads
            );
            
            filtered.sort((a, b) => {
                if (sortBy === 'callers') return b.callers.length - a.callers.length;
                if (sortBy === 'reads') return b.totalReads - a.totalReads;
                if (sortBy === 'name') return a.fusionName.localeCompare(b.fusionName);
                return 0;
            });
            
            return filtered;
        }
        
        // Render fusion table
        function renderFusionTable() {
            const fusions = getFilteredFusions();
            const allCallers = ['Arriba', 'Squid', 'Pizzly', 'FusionCatcher', 'STAR-Fusion'];
            
            let html = '<div class="table-container"><table>';
            html += '<thead><tr>';
            html += '<th>Fusion</th>';
            html += '<th class="text-center">Callers</th>';
            allCallers.forEach(caller => {
                html += `<th colspan="2" class="text-center">${caller}</th>`;
            });
            html += '<th class="text-right">Total Reads</th>';
            html += '</tr><tr>';
            html += '<th></th><th></th>';
            allCallers.forEach(() => {
                html += '<th class="text-center">Split</th>';
                html += '<th class="text-center">Spanning</th>';
            });
            html += '<th></th>';
            html += '</tr></thead>';
            html += '<tbody>';
            
            fusions.forEach((fusion, idx) => {
                html += `<tr class="fusion-row" onclick="showFusionDetail(${idx})">`;
                html += `<td class="font-bold text-primary">${fusion.fusionName}</td>`;
                html += `<td class="text-center"><span class="badge">${fusion.callers.length}</span></td>`;
                
                allCallers.forEach(caller => {
                    const bp = fusion.breakpoints.find(b => b.caller === caller);
                    const split = bp ? (bp.splitReads || '-') : '-';
                    const spanning = bp ? (bp.discordantMates || '-') : '-';
                    html += `<td class="text-center">${split}</td>`;
                    html += `<td class="text-center">${spanning}</td>`;
                });
                
                html += `<td class="text-right font-bold">${fusion.totalReads.toLocaleString()}</td>`;
                html += '</tr>';
            });
            
            html += '</tbody></table></div>';
            document.getElementById('fusionsTable').innerHTML = html;
        }
        
        // Render coverage table
        function renderCoverageTable() {
            let html = '<div class="table-container"><table>';
            html += '<thead><tr>';
            html += '<th>Chromosome</th>';
            html += '<th>Start</th>';
            html += '<th>End</th>';
            html += '<th>Target</th>';
            html += '<th class="text-center">Strand</th>';
            html += '<th class="text-right">Coverage</th>';
            html += '</tr></thead>';
            html += '<tbody>';
            
            coverageData.forEach(row => {
                html += '<tr>';
                html += `<td>${row.chr}</td>`;
                html += `<td>${row.start.toLocaleString()}</td>`;
                html += `<td>${row.end.toLocaleString()}</td>`;
                html += `<td>${row.target}</td>`;
                html += `<td class="text-center">${row.strand}</td>`;
                html += `<td class="text-right font-bold">${row.coverage.toLocaleString()}</td>`;
                html += '</tr>';
            });
            
            html += '</tbody></table></div>';
            document.getElementById('coverageTable').innerHTML = html;
        }
        
        // Show fusion detail inline
        function showFusionDetail(index) {
            const fusions = getFilteredFusions();
            currentFusion = fusions[index];
            
            document.getElementById('detailTitle').innerHTML = 
                `Fusion Details: <span class="text-primary">${currentFusion.fusionName}</span>`;
            
            document.getElementById('detailStats').innerHTML = `
                <div class="stat-row">
                    <div class="stat-card">
                        <div class="stat-value">${currentFusion.callers.length}</div>
                        <div class="stat-label">Callers</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">${currentFusion.totalReads.toLocaleString()}</div>
                        <div class="stat-label">Total Reads</div>
                    </div>
                </div>
            `;

            
            document.getElementById('callersList').textContent = 
                `Detected by: ${currentFusion.callers.join(', ')}`;
            
            renderBreakpointsTable();
            renderCircosPlot();
            
            const detailElement = document.getElementById('fusionDetail');
            detailElement.classList.add('active');
            
            // Auto-scroll to fusion detail
            setTimeout(() => {
                detailElement.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }, 100);
        }
        
        // Render breakpoints table
        function renderBreakpointsTable() {
            let html = '<div class="table-container"><table>';
            html += '<thead><tr>';
            html += '<th>Caller</th>';
            html += '<th>Breakpoint 1</th>';
            html += '<th>Breakpoint 2</th>';
            html += '<th class="text-center">Split Reads</th>';
            html += '<th class="text-center">Discordant Mates</th>';
            html += '<th style="max-width: 300px;">Annotation</th>';
            html += '<th class="text-center">Confidence</th>';
            html += '<th>Coding Frames</th>';
            html += '<th>Type</th>';
            html += '<th>Retained Protein Domains</th>';
            html += '</tr></thead>';
            html += '<tbody>';
            
            currentFusion.breakpoints.forEach(bp => {
                html += '<tr>';
                html += `<td><span class="badge secondary">${bp.caller}</span></td>`;
                html += `<td style="font-family: monospace; font-size: 0.875rem">${bp.breakpoint1}</td>`;
                html += `<td style="font-family: monospace; font-size: 0.875rem">${bp.breakpoint2}</td>`;
                html += `<td class="text-center">${bp.splitReads || '-'}</td>`;
                html += `<td class="text-center">${bp.discordantMates || '-'}</td>`;
                html += `<td style="max-width: 300px; word-wrap: break-word; white-space: normal;">${bp.annotation || '-'}</td>`;
                html += `<td class="text-center">${bp.confidence ? `<span class="badge">${bp.confidence}</span>` : '-'}</td>`;
                html += `<td>${bp.codingFrames || '-'}</td>`;
                html += `<td>${bp.type || '-'}</td>`;
                html += `<td>${bp.retainedProteinDomains || '-'}</td>`;
                html += '</tr>';
            });
            
            html += '</tbody></table></div>';
            document.getElementById('breakpointsTable').innerHTML = html;
        }
        
        // Render circos plot
        function renderCircosPlot() {
            const container = document.getElementById('circosPlot');
            container.innerHTML = '';
            
            const width = 700;
            const height = 700;
            const radius = Math.min(width, height) / 2 - 70;
            
            const svg = d3.select('#circosPlot')
                .append('svg')
                .attr('width', width)
                .attr('height', height);
            
            const g = svg.append('g')
                .attr('transform', `translate(${width/2},${height/2})`);
            
            // Group cytobands by chromosome
            const chrData = new Map();
            cytobandData.forEach(band => {
                if (!chrData.has(band.chr)) {
                    chrData.set(band.chr, { bands: [], length: 0 });
                }
                const data = chrData.get(band.chr);
                data.bands.push(band);
                data.length = Math.max(data.length, band.end);
            });
            
            // Sort chromosomes
            const chrArray = Array.from(chrData.keys()).sort((a, b) => {
                const aNum = a.replace('chr', '');
                const bNum = b.replace('chr', '');
                if (aNum === 'X') return 1;
                if (bNum === 'X') return -1;
                if (aNum === 'Y') return 1;
                if (bNum === 'Y') return -1;
                return parseInt(aNum) - parseInt(bNum);
            });
            
            const totalGenomeLength = Array.from(chrData.values())
                .reduce((sum, chr) => sum + chr.length, 0);
            
            const stainColors = {
                'gneg': '#FFFFFF',
                'gpos25': '#C0C0C0',
                'gpos50': '#808080',
                'gpos75': '#404040',
                'gpos100': '#000000',
                'acen': '#8B4513',
                'gvar': '#4682B4'
            };
            
            let currentAngle = 0;
            const chrPositions = new Map();
            
            // Chromosome colors for inner circle boundary identification
            const chrColors = ['#A8D5E2', '#C8E6C9', '#FFE0B2', '#F8BBD0', '#D1C4E9', '#FFCCBC'];
            
            // Draw inner circle with chromosome color segments
            chrArray.forEach((chr, chrIdx) => {
                const chrInfo = chrData.get(chr);
                const chrAngleSpan = (chrInfo.length / totalGenomeLength) * 2 * Math.PI;
                const startAngle = currentAngle;
                const endAngle = currentAngle + chrAngleSpan;
                
                chrPositions.set(chr, { startAngle, endAngle, length: chrInfo.length });
                
                // Draw inner colored arc for chromosome identification
                const innerArc = d3.arc()
                    .innerRadius(radius - 50)
                    .outerRadius(radius - 35)
                    .startAngle(startAngle)
                    .endAngle(endAngle);
                
                g.append('path')
                    .attr('d', innerArc)
                    .attr('fill', chrColors[chrIdx % chrColors.length])
                    .attr('stroke', '#666')
                    .attr('stroke-width', 1.5);
                
                // Draw chromosome underlay for cytoband background
                const underlayArc = d3.arc()
                    .innerRadius(radius - 30)
                    .outerRadius(radius)
                    .startAngle(startAngle)
                    .endAngle(endAngle);
                
                g.append('path')
                    .attr('d', underlayArc)
                    .attr('fill', '#f5f5f5')
                    .attr('stroke', '#999')
                    .attr('stroke-width', 1.5)
                    .attr('opacity', 0.3);
                
                // Draw cytobands
                chrInfo.bands.forEach(band => {
                    const bandStartAngle = startAngle + (band.start / chrInfo.length) * chrAngleSpan;
                    const bandEndAngle = startAngle + (band.end / chrInfo.length) * chrAngleSpan;
                    
                    const arc = d3.arc()
                        .innerRadius(radius - 30)
                        .outerRadius(radius)
                        .startAngle(bandStartAngle)
                        .endAngle(bandEndAngle);
                    
                    g.append('path')
                        .attr('d', arc)
                        .attr('fill', stainColors[band.stain] || '#CCCCCC')
                        .attr('stroke', '#ddd')
                        .attr('stroke-width', 0.3);
                });
                
                // Add chromosome label
                const labelAngle = (startAngle + endAngle) / 2;
                const labelRadius = radius + 25;
                const x = Math.cos(labelAngle) * labelRadius;
                const y = Math.sin(labelAngle) * labelRadius;
                
                g.append('text')
		    .attr('x', x)
		    .attr('y', y)
		    .attr('text-anchor', 'middle')
		    .attr('dominant-baseline', 'middle')
		    .attr('fill', '#333')
		    .attr('font-size', '11px')
		    .attr('font-weight', 'bold')
		    .text(chr.replace('chr', ''));
				
                currentAngle = endAngle;
            });
            
            // Draw fusion links
            currentFusion.breakpoints.forEach((bp, idx) => {
                // Normalize chromosome names so “1” → “chr1”
		const chr1 = bp.breakpoint1.split(':')[0].startsWith('chr')
		    ? bp.breakpoint1.split(':')[0]
		    : 'chr' + bp.breakpoint1.split(':')[0];

		const chr2 = bp.breakpoint2.split(':')[0].startsWith('chr')
		    ? bp.breakpoint2.split(':')[0]
		    : 'chr' + bp.breakpoint2.split(':')[0];

                const pos1 = Number(bp.breakpoint1.split(':')[1]);
		const pos2 = Number(bp.breakpoint2.split(':')[1]);

                
                const chrInfo1 = chrPositions.get(chr1);
                const chrInfo2 = chrPositions.get(chr2);
                
                if (!chrInfo1 || !chrInfo2 || isNaN(pos1) || isNaN(pos2)) return;
                
                const angle1 = chrInfo1.startAngle + 
                    (pos1 / chrInfo1.length) * (chrInfo1.endAngle - chrInfo1.startAngle);
                const angle2 = chrInfo2.startAngle + 
                    (pos2 / chrInfo2.length) * (chrInfo2.endAngle - chrInfo2.startAngle);
                
                const linkRadius = radius - 40;
                const x1 = Math.cos(angle1) * linkRadius;
                const y1 = Math.sin(angle1) * linkRadius;
                const x2 = Math.cos(angle2) * linkRadius;
                const y2 = Math.sin(angle2) * linkRadius;
                
                const path = d3.path();
                path.moveTo(x1, y1);
                path.quadraticCurveTo(0, 0, x2, y2);
                
                g.append('path')
                    .attr('d', path.toString())
                    .attr('fill', 'none')
                    .attr('stroke', `hsl(${(idx * 60) % 360}, 85%, 55%)`)
                    .attr('stroke-width', 2.5)
                    .attr('opacity', 0.7);
                
                g.append('circle')
                    .attr('cx', x1)
                    .attr('cy', y1)
                    .attr('r', 4)
                    .attr('fill', 'hsl(217, 91%, 60%)')
                    .attr('stroke', '#fff')
                    .attr('stroke-width', 1);
                
                g.append('circle')
                    .attr('cx', x2)
                    .attr('cy', y2)
                    .attr('r', 4)
                    .attr('fill', 'hsl(0, 84%, 60%)')
                    .attr('stroke', '#fff')
                    .attr('stroke-width', 1);
            });
            
            // Add gene labels
            g.append('text')
                .attr('text-anchor', 'middle')
                .attr('y', -10)
                .attr('fill', '#333')
                .attr('font-size', '18px')
                .attr('font-weight', 'bold')
                .text(currentFusion.gene1);
            
            g.append('text')
                .attr('text-anchor', 'middle')
                .attr('y', 15)
                .attr('fill', '#333')
                .attr('font-size', '14px')
                .text('--');
            
            g.append('text')
                .attr('text-anchor', 'middle')
                .attr('y', 40)
                .attr('fill', '#333')
                .attr('font-size', '18px')
                .attr('font-weight', 'bold')
                .text(currentFusion.gene2);
        }
        
        // Event listeners
        document.getElementById('minCallers').addEventListener('change', renderFusionTable);
        document.getElementById('minReads').addEventListener('input', renderFusionTable);
        document.getElementById('sortBy').addEventListener('change', renderFusionTable);
        
        // Initial render
        updateSummaryCards();
        renderFusionTable();
        renderCoverageTable();
    </script>
</body>
</html>'''
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_template)

def main():
    parser = argparse.ArgumentParser(
        description="Generate gene fusion analysis dashboard HTML report."
    )

    parser.add_argument(
        "--fusions", "-f",
        required=True,
        help="Excel file containing fusion caller outputs (.xlsx)"
    )

    parser.add_argument(
        "--cytoband", "-c",
        required=True,
        help="Cytoband annotation file (hg38)"
    )

    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Optional: output HTML filename. If omitted, auto-generated as <sample>_dashboard.html"
    )

    args = parser.parse_args()

    excel_file = args.fusions
    cytoband_file = args.cytoband

    # Auto-generate output filename if not provided
    if args.output:
        output_file = args.output
    else:
        sample_name = Path(excel_file).stem
        output_file = f"{sample_name}_dashboard.html"

    sample_name = Path(excel_file).stem  # still used inside HTML

    # Validate input files
    if not Path(excel_file).exists():
        print(f"Error: Excel file not found: {excel_file}")
        sys.exit(1)

    if not Path(cytoband_file).exists():
        print(f"Error: Cytoband file not found: {cytoband_file}")
        sys.exit(1)

    print("Parsing coverage data...")
    coverage_data = parse_coverage_sheet(excel_file)
    print(f"Found {len(coverage_data)} coverage records")

    print("\nLoading fusion caller sheets...")

    fusion_dfs = {}

    # Arriba
    try:
        fusion_dfs["Arriba"] = pd.read_excel(excel_file, sheet_name=".arriba.fusions")
    except:
        print("Warning: .arriba.fusions sheet missing")

    # Squid
    try:
        fusion_dfs["Squid"] = pd.read_excel(excel_file, sheet_name=".squid.fusions.annotated")
    except:
        print("Warning: .squid.fusions.annotated sheet missing")

    # Pizzly
    try:
        fusion_dfs["Pizzly"] = pd.read_excel(excel_file, sheet_name=".pizzly")
    except:
        print("Warning: .pizzly sheet missing")

    # FusionCatcher
    try:
        fusion_dfs["FusionCatcher"] = pd.read_excel(excel_file, sheet_name=".fusioncatcher.fusion-genes")
    except:
        print("Warning: .fusioncatcher.fusion-genes sheet missing")

    # STAR-Fusion
    try:
        fusion_dfs["STAR-Fusion"] = pd.read_excel(excel_file, sheet_name=".starfusion.fusion_predictions")
    except:
        print("Warning: .starfusion.fusion_predictions sheet missing")

    print("\nLoaded Fusion Sheets:")
    for k, v in fusion_dfs.items():
        print(f"{k}: {len(v)} rows")

    print("\nParsing fusion data...")
    fusion_data = parse_fusion_callers(fusion_dfs)
    print(f"Found {len(fusion_data)} unique fusions")

    print("\nLoading cytoband data...")
    cytoband_data = load_cytoband_data(cytoband_file)
    print(f"Loaded {len(cytoband_data)} cytoband records")

    print("\nGenerating HTML dashboard...")
    generate_html(coverage_data, fusion_data, cytoband_data, output_file, sample_name)
    print(f"Dashboard generated successfully: {output_file}")


if __name__ == '__main__':
    main()
