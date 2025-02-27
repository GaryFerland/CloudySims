import ads
import time
import matplotlib.pyplot as plt
from collections import defaultdict

# Set ADS API token (Replace with your own ADS token)
ads.config.token = "YOUR_ADS_API_TOKEN`"

# **Step 1: Define Bibcodes for Cloudy Versions**
bibcode_queries = {
    "Cloudy C90": "1998PASP..110..761F",
    "Cloudy C13": "2013RMxAA..49..137F",
    "Cloudy C17": "2017RMxAA..53..385F",
    "Cloudy C23.01": "2023RNAAS...7..246G",
    "Cloudy C23": "2023RMxAA..59..327C"
}

# **Step 2: Define Hazy Bibcodes**
hazy_bibcodes = [
    "1988hbic.book.....F", "2006hbic.book.....F", "1987hbic.book.....F",
    "2000hbic.book.....F", "2004hbic.book.....F", "1989hbic.book.....F",
    "2002hbic.book.....F", "1990hbic.book.....F", "1991hbic.book.....F",
    "1993hbic.book.....F", "1996hbic.book.....F", "1983hbic.book.....F", "2005hbic.book.....F"
]

# **Step 3: Fetch Citations from ADS**
def fetch_unique_citations(query, refereed_only=False):
    """ Fetch citations from ADS with a fallback for refereed papers. """
    citing_papers = set()
    start = 0
    rows = 2000

    while True:
        try:
            q = query + " property:refereed" if refereed_only else query
            papers = ads.SearchQuery(q=q, fl=["year", "bibcode"], rows=rows, start=start)
            papers_list = list(papers)

            # If no refereed results, fetch non-refereed citations
            if not papers_list and refereed_only:
                papers = ads.SearchQuery(q=query, fl=["year", "bibcode"], rows=rows, start=start)
                papers_list = list(papers)

            if not papers_list:
                break

            for paper in papers_list:
                if paper.year and 1980 <= int(paper.year) <= 2025:
                    citing_papers.add((paper.bibcode, int(paper.year)))

            start += rows

        except Exception as e:
            print(f"Error fetching citations: {e}")
            break

    return citing_papers

# **Step 4: Fetch Hazy Citations Using ADS Metrics API**
def fetch_hazy_citations(bibcodes):
    """ Fetch citations for Hazy books using ADS Metrics API. """
    hazy_citations_by_year = defaultdict(int)

    for bibcode in bibcodes:
        try:
            metrics = ads.MetricsQuery(bibcodes=[bibcode]).execute()
            histograms = metrics.get("histograms", {}).get("citations", {})

            for category in ["refereed to refereed", "refereed to nonrefereed", "nonrefereed to refereed", "nonrefereed to nonrefereed"]:
                if category in histograms:
                    for year_str, count in histograms[category].items():
                        try:
                            year = int(year_str)
                            if 1980 <= year <= 2025:
                                hazy_citations_by_year[year] += count
                        except ValueError:
                            pass
        except Exception as e:
            print(f"Error fetching Hazy citations: {e}")
            continue

    return hazy_citations_by_year

# **Step 5: Store Cloudy Citation Data**
data = {}
for version, bibcode in bibcode_queries.items():
    papers = fetch_unique_citations(f"citations({bibcode})")
    citations_by_year = defaultdict(int)
    
    for _, year in papers:
        citations_by_year[year] += 1
    
    if version in ["Cloudy C23", "Cloudy C23.01"]:
        version = "Cloudy C23"

    if version not in data:
        data[version] = defaultdict(int)

    for year, count in citations_by_year.items():
        data[version][year] += count

# **Step 6: Fetch Hazy Citations and Store in Data**
hazy_data = fetch_hazy_citations(hazy_bibcodes)
data["Hazy"] = hazy_data

# **Step 7: Fetch Refereed Citations for Cloudy Versions**
refereed_citations_by_year = defaultdict(int)
for version, bibcode in bibcode_queries.items():
    papers = fetch_unique_citations(f"citations({bibcode})", refereed_only=True)
    for _, year in papers:
        refereed_citations_by_year[year] += 1

# **Step 8: Ensure Refereed Citations = Hazy + Cloudy Refereed**
for year in sorted(hazy_data):
    refereed_citations_by_year[year] += hazy_data[year]  # Add Hazy citations to refereed count

# **Step 9: Prepare Data for Plotting**
years = sorted(set(year for version_data in data.values() for year in version_data))
ordered_versions = ["Hazy", "Cloudy C90", "Cloudy C13", "Cloudy C17", "Cloudy C23"]
version_colors = ["blue", "orange", "green", "red", "purple"]

# **Step 10: Plot Stacked Bar Chart**
plt.figure(figsize=(12, 8))
bottom = [0] * len(years)

for version, color in zip(ordered_versions, version_colors):
    counts = [data[version].get(year, 0) for year in years]
    plt.bar(years, counts, bottom=bottom, label=version, color=color)
    bottom = [b + c for b, c in zip(bottom, counts)]

# **Step 11: Overlay Refereed Citations as a Line Plot**
refereed_counts = [refereed_citations_by_year.get(year, 0) for year in years]
plt.plot(years, refereed_counts, label="Refereed Papers", color="black", marker="o", linewidth=2)

# **Step 12: Add Labels, Legend, and Title**
plt.xlabel("Year")
plt.ylabel("Citations")
plt.title("Citations of Cloudy Versions and Hazy Over Time (Refereed Papers Included)")
plt.grid(which="both", linestyle="--", linewidth=0.5)
plt.xticks(ticks=years, labels=years, rotation=90)
plt.legend()
plt.tight_layout()
plt.show()

# **Step 13: Print Final Summary**
print("\n📊 **Total Citations by Version:**")
for version in ordered_versions:
    total_citations = sum(data[version].values())
    print(f"{version}: {total_citations} citations")

print("\n📊 **Total Refereed Citations Per Year:**")
for year in sorted(refereed_citations_by_year):
    print(f"{year}: {refereed_citations_by_year[year]} refereed citations")

