

import time
import os
import csv
from io import StringIO

from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import pandas as pd


class NCBIRetriever:
    """
    Handles interaction with the NCBI GenBank database.
    """

    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        self.webenv = None
        self.query_key = None
        self.count = 0
        self.organism_name = ""

        # Set up Entrez
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        Entrez.tool = 'BioScriptEx10_sXXXXX'  # Please change sXXXXX

    def search_taxid(self, taxid, min_len=None, max_len=None):
        """
        Search for all records associated with a taxonomic ID,
        optionally filtered by sequence length.
        """
        print(f"Searching for records with taxID: {taxid}")
        if min_len is not None or max_len is not None:
            len_filter_str = " and "
            if min_len is not None and max_len is not None:
                len_filter_str += f"{min_len}:{max_len}[SLEN]"
            elif min_len is not None:
                len_filter_str += f"{min_len}:99999999999[SLEN]"  # Effectively min_len to infinity
            else:  # only max_len is not None
                len_filter_str += f"1:{max_len}[SLEN]"  # Effectively 1 to max_len
            print(f"Applying sequence length filter: {len_filter_str.replace('[SLEN]', '').replace(' and ', '')}")
        else:
            len_filter_str = ""

        try:
            # Get taxonomy information first
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            if not records:
                print(f"Could not retrieve taxonomy information for TaxID {taxid}.")
                return None
            self.organism_name = records[0]["ScientificName"]
            print(f"Organism: {self.organism_name} (TaxID: {taxid})")

            # Search for records
            search_term = f"txid{taxid}[Organism]{len_filter_str}"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", idtype="acc")
            search_results = Entrez.read(handle)
            handle.close()

            self.count = int(search_results["Count"])

            if self.count == 0:
                print(f"No records found for {self.organism_name} with the specified criteria.")
                return 0

            print(f"Found {self.count} records matching criteria.")

            # Store search results for later use
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]

            return self.count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_all_sequences_data(self, max_records_to_process=None):
        """
        Fetch all records using the stored search results, parse them,
        and extract accession, length, and description.
        """
        if not self.webenv or not self.query_key:
            print("No search results to fetch. Run search_taxid() first.")
            return []

        all_sequence_details = []

        records_to_retrieve = self.count
        if max_records_to_process is not None and max_records_to_process < self.count:
            print(f"Will process a maximum of {max_records_to_process} records out of {self.count} found.")
            records_to_retrieve = max_records_to_process

        batch_size = 200  # NCBI recommended batch size
        retrieved_count = 0

        for start_index in range(0, records_to_retrieve, batch_size):
            current_batch_size = min(batch_size, records_to_retrieve - start_index)
            try:
                print(f"Fetching records {start_index + 1} to {start_index + current_batch_size}...")
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start_index,
                    retmax=current_batch_size,
                    webenv=self.webenv,
                    query_key=self.query_key,
                    idtype="acc"
                )
                records_text = handle.read()
                handle.close()

                # Parse the GenBank records
                # Using StringIO to treat the string as a file
                gb_iter = SeqIO.parse(StringIO(records_text), "gb")
                for record in gb_iter:
                    all_sequence_details.append({
                        'accession': record.id,
                        'length': len(record.seq),
                        'description': record.description
                    })
                    retrieved_count += 1

                # Respect NCBI rate limits (10 req/sec with API key, 3 req/sec without)
                time.sleep(0.11 if self.api_key else 0.34)

            except Exception as e:
                print(f"Error fetching or parsing records batch (start: {start_index}): {e}")
                # Optionally, decide if you want to stop or continue with next batch
                time.sleep(1)  # Wait a bit longer after an error

        print(f"Successfully fetched and parsed {retrieved_count} records.")
        return all_sequence_details


def generate_csv_report(data, taxid, organism_name, min_len, max_len):
    """
    Generates a CSV file containing GenBank accession, sequence length, and description.
    """
    if not data:
        print("No data to generate CSV report.")
        return

    # Create a filename that reflects the query parameters
    len_suffix = ""
    if min_len is not None or max_len is not None:
        min_str = f"_min{min_len}" if min_len is not None else "_minAny"
        max_str = f"_max{max_len}" if max_len is not None else "_maxAny"
        len_suffix = f"{min_str}{max_str}"

    filename = f"taxid_{taxid}_{organism_name.replace(' ', '_')}{len_suffix}_report.csv"

    try:
        with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['accession', 'length', 'description']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        print(f"CSV report generated: {filename}")
    except IOError as e:
        print(f"Error writing CSV file {filename}: {e}")


def generate_length_chart(data, taxid, organism_name, min_len, max_len):
    """
    Generates a line chart of sequence length vs. GenBank accession, sorted by length.
    Saves the chart as a PNG file.
    """
    if not data:
        print("No data to generate chart.")
        return

    df = pd.DataFrame(data)
    if df.empty:
        print("DataFrame is empty, cannot generate chart.")
        return

    df_sorted = df.sort_values(by='length', ascending=False)

    plt.figure(figsize=(12, 7))  # Adjust figure size as needed

    # Using a range for x-ticks if accessions are too numerous or too long
    if len(df_sorted['accession']) > 50:  # Heuristic for too many labels
        x_values = range(len(df_sorted['accession']))
        plt.plot(x_values, df_sorted['length'])
        plt.xticks(x_values[::max(1, len(x_values) // 20)], df_sorted['accession'][::max(1, len(x_values) // 20)],
                   rotation=90, fontsize=8)  # Show limited ticks
        plt.xlabel("GenBank Accession (sampled)")
    else:
        plt.plot(df_sorted['accession'], df_sorted['length'])
        plt.xlabel("GenBank Accession Number")
        plt.xticks(rotation=90, fontsize=8)  # Rotate for readability

    plt.ylabel("Sequence Length (bp)")
    title_org_name = organism_name if len(organism_name) < 40 else organism_name[:37] + "..."
    plt.title(f"Sequence Length Distribution for {title_org_name} (TaxID: {taxid})")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()  # Adjust layout to prevent labels from being cut off

    # Create a filename that reflects the query parameters
    len_suffix = ""
    if min_len is not None or max_len is not None:
        min_str = f"_min{min_len}" if min_len is not None else "_minAny"
        max_str = f"_max{max_len}" if max_len is not None else "_maxAny"
        len_suffix = f"{min_str}{max_str}"

    filename = f"taxid_{taxid}_{organism_name.replace(' ', '_')}{len_suffix}_length_chart.png"

    try:
        plt.savefig(filename)
        print(f"Chart saved as: {filename}")
    except IOError as e:
        print(f"Error saving chart file {filename}: {e}")
    finally:
        plt.close()  # Close the plot to free memory


def get_length_input(prompt_message):
    """Helper function to get and validate length input."""
    while True:
        try:
            val_str = input(prompt_message).strip()
            if not val_str:  # User pressed Enter for no limit
                return None
            val = int(val_str)
            if val < 0:
                print("Length cannot be negative. Please enter a valid number or leave blank.")
                continue
            return val
        except ValueError:
            print("Invalid input. Please enter a whole number or leave blank.")


def main():
    """
    Main function to drive the script.
    """
    print("NCBI GenBank Data Retriever and Analyzer")
    print("-" * 40)

    email = input("Enter your email address for NCBI: ").strip()
    api_key = input("Enter your NCBI API key (optional, press Enter to skip): ").strip()
    if not api_key:
        api_key = None
        print("Proceeding without API key. Rate limits will be 3 requests/second.")
    else:
        print("Using API key. Rate limits will be 10 requests/second.")

    retriever = NCBIRetriever(email, api_key)

    taxid = input("Enter taxonomic ID (taxid) of the organism (e.g., 9606 for Homo sapiens): ").strip()
    if not taxid.isdigit():
        print("Invalid TaxID. Please enter a numerical ID.")
        return

    print("\nSpecify sequence length thresholds (press Enter for no limit):")
    min_len = get_length_input("Enter minimum sequence length (e.g., 1000): ")
    max_len = get_length_input("Enter maximum sequence length (e.g., 50000): ")

    if min_len is not None and max_len is not None and min_len > max_len:
        print("Error: Minimum length cannot be greater than maximum length.")
        return

    # Search for records
    count = retriever.search_taxid(taxid, min_len, max_len)

    if count is None:  # Indicates an error during search
        print("Search operation failed. Exiting.")
        return
    if count == 0:
        print("No records found matching your criteria. Exiting.")
        return
    process_limit = count
    if count > 500:  # Arbitrary threshold to ask for confirmation if many records
        try:
            limit_input = input(
                f"Found {count} records. Enter maximum number to process (or press Enter for all {count}): ").strip()
            if limit_input:
                process_limit = int(limit_input)
                if process_limit <= 0:
                    print("Number of records to process must be positive. Defaulting to all.")
                    process_limit = count
                elif process_limit > count:
                    print(f"Input ({process_limit}) exceeds found records ({count}). Processing all {count} records.")
                    process_limit = count
                else:
                    print(f"Will process up to {process_limit} records.")
            else:
                print(f"Will process all {count} records.")
        except ValueError:
            print("Invalid input for record limit. Defaulting to all.")
            process_limit = count

    # Fetch sequence details
    print(f"\nFetching details for up to {process_limit} records...")
    sequence_details = retriever.fetch_all_sequences_data(max_records_to_process=process_limit)

    if not sequence_details:
        print("No sequence details could be retrieved. Exiting.")
        return

    # Generate CSV report
    print("\nGenerating CSV report...")
    generate_csv_report(sequence_details, taxid, retriever.organism_name, min_len, max_len)

    # Generate Data Visualization
    print("\nGenerating data visualization chart...")
    generate_length_chart(sequence_details, taxid, retriever.organism_name, min_len, max_len)

    print("\nScript finished.")

if __name__ == "__main__":
    main()