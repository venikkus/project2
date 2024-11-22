import re
import pandas as pd


def process_vcf(vcf_file):
    """Parse VCF file to extract position, allele and frequency data."""
    data = {
        "Position": [],
        "Ref": [],
        "Alt": [],
        "Freq": []
    }

    with open(vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            ref, alt, pos = columns[3], columns[4], columns[1]

            match = re.search(r"([\d,\.]+%)", columns[9])
            freq = float(match.group(1).strip('%').replace(',', '.')) \
                if match else None

            data["Position"].append(pos)
            data["Ref"].append(ref)
            data["Alt"].append(alt)
            data["Freq"].append(freq)

    return pd.DataFrame(data)


def calculate_thresholds(df_reference) -> float:
    """Bounds calsulations based on reference file."""
    filtered_ref = df_reference[df_reference['Freq'] <= 0.9]
    mean_freq = filtered_ref['Freq'].mean()
    std_freq = filtered_ref['Freq'].std()
    return mean_freq, std_freq, mean_freq + 3 * std_freq


def vcf_to_table(vcf_file1, vcf_file2) -> float:
    """VCF-file filtered based on bounds."""
    df_roommate = process_vcf(vcf_file1)
    df_reference = process_vcf(vcf_file2)

    mean_freq_ref, std_freq_ref, threshold = \
        calculate_thresholds(df_reference)

    filtered_df = df_roommate[
        (df_roommate['Freq'] > threshold) & (df_roommate['Freq'] < 90)
        ]
    return filtered_df, mean_freq_ref, std_freq_ref


filenames = []
print("Write filename of experiment sample and then filenames for control samples "
      "(empty string to complete): ")

while True:
    name = input()
    if name == "":
        break
    filenames.append(name)

filename = filenames[0]
reference_files = filenames[1:]

all_filtered_values = []
means_ref, stds_ref = [], []

for ref_file in reference_files:
    filtered_values, mean_ref, std_ref = vcf_to_table(filename, ref_file)
    all_filtered_values.append(filtered_values)
    means_ref.append(mean_ref)
    stds_ref.append(std_ref)

combined_filtered = pd.concat(all_filtered_values, ignore_index=True)

mean_std = sum(stds_ref) / len(stds_ref)
mean_mean = sum(means_ref) / len(means_ref)

filtered_above_mean_std = combined_filtered[
    combined_filtered['Freq'] > 3 * mean_std
].drop_duplicates()

print("\nFiltered values that exceed the mean of std deviations:")
print(filtered_above_mean_std)

print("\nSummary of Means and Standard Deviations:")
for i, ref_file in enumerate(reference_files):
    print(f"{ref_file}: Mean = {means_ref[i]:.4f}, Std = {stds_ref[i]:.4f}")

print(f"\nMean = {mean_mean:.4f}, Std = {mean_std:.4f} \
      \nMean + Std * 3 = {mean_mean + (mean_std * 3):.4f}")
