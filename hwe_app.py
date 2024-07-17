import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, fisher_exact, chi2_contingency
from statsmodels.stats.power import GofChisquarePower
from statsmodels.stats.contingency_tables import Table
from sklearn.linear_model import LogisticRegression
import pandas as pd
from itertools import combinations
import requests

# Title and introduction
st.title("Hardy-Weinberg Equilibrium Analysis App for Clinical Genetics")

# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Advanced Statistical Tests", "Allele Frequency Evolution", "Statistical Power Analysis", "Clinical Interpretation"])

# Home Page
if page == "Home":
    st.header("Welcome to the HWE Analysis App")
    st.write("Use the navigation on the left to access different features.")
    st.write("""
    This app allows you to perform Hardy-Weinberg Equilibrium (HWE) analysis with advanced features including:
    - Allele Frequency Evolution Simulation
    - Statistical Power Analysis
    - Advanced Statistical Tests
    - Clinical Interpretation of HWE Deviations
    """)
    
# Function to calculate expected counts under HWE
def hwe_expected(obs_hom_var, obs_het, obs_hom_ref):
    total = obs_hom_var + obs_het + obs_hom_ref
    total_alleles = 2 * total
    
    num_a_alleles = 2 * obs_hom_var + obs_het
    num_A_alleles = 2 * obs_hom_ref + obs_het
    
    q = num_a_alleles / total_alleles
    p = num_A_alleles / total_alleles
    
    exp_hom_var = (q**2) * total
    exp_het = 2 * p * q * total
    exp_hom_ref = (p**2) * total
    return exp_hom_var, exp_het, exp_hom_ref

# Function for allele frequency simulation
def simulate_allele_freq(p0, selection_coef, migration_rate, population_size, generations):
    freqs = [p0]
    p = p0
    for _ in range(generations):
        q = 1 - p
        p_next = p * (1 - selection_coef) + migration_rate * q
        p_next /= (p_next + q * (1 - selection_coef))
        freqs.append(p_next)
        p = p_next
    return freqs

# Function for power calculation
def calculate_power(effect_size, nobs, alpha):
    power_analysis = GofChisquarePower()
    power = power_analysis.solve_power(effect_size=effect_size, nobs=nobs, alpha=alpha)
    return power

# Function for exact test
def exact_test(obs):
    _, p_value = fisher_exact(obs)
    return p_value

# Function for chi-square test with multiple alleles
def chi_square_test(obs):
    chi2, p_value, _, _ = chi2_contingency(obs)
    return chi2, p_value

# Function to fetch gene and variant information
def fetch_variant_info(variant):
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{variant}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return None

# Monte Carlo simulation for Fisher's exact test on a 2x3 table
def monte_carlo_fishers_exact(table, num_simulations=10000):
    obs = table.sum(axis=0)
    n = table.sum()
    count = 0
    for _ in range(num_simulations):
        shuffled = np.random.choice(obs, size=n, replace=True)
        simulated_table = np.array([
            [shuffled[0], shuffled[1], shuffled[2]],
            [shuffled[3], shuffled[4], shuffled[5]]
        ])
        if fisher_exact(simulated_table, alternative='two-sided')[1] >= fisher_exact(table, alternative='two-sided')[1]:
            count += 1
    return count / num_simulations

# Permutation test for a 2x3 table
def permutation_test(table, num_permutations=10000):
    observed_stat, _, _, _ = chi2_contingency(table)
    count = 0
    combined = table.flatten()
    for _ in range(num_permutations):
        np.random.shuffle(combined)
        permuted_table = combined.reshape(table.shape)
        permuted_stat, _, _, _ = chi2_contingency(permuted_table)
        if permuted_stat >= observed_stat:
            count += 1
    return count / num_permutations

if page == "Advanced Statistical Tests":
    st.header("Advanced Statistical Tests")
    
    # Explanation of HWE Expected Calculation
    st.write("""
    ### How Expected Genotype Frequencies are Calculated
    Under Hardy-Weinberg Equilibrium (HWE), the expected genotype frequencies are calculated using the allele frequencies.
    If p is the frequency of the major allele and q is the frequency of the minor allele (where p + q = 1), then:
    - Expected frequency of homozygous major (AA) = p^2
    - Expected frequency of heterozygous (Aa) = 2pq
    - Expected frequency of homozygous minor (aa) = q^2
    
    These frequencies are then multiplied by the total number of individuals to get the expected genotype counts.
    """)

    # User inputs
    st.write("### Enter the observed genotype counts:")
    homozygous_variant = st.number_input("Homozygous Variant", 1, 100000, 1, key="hom_var")
    heterozygous = st.number_input("Heterozygous", 1, 100000, 1, key="het")
    homozygous_reference = st.number_input("Homozygous Reference", 1, 100000, 1, key="hom_ref")

    # Calculate expected counts
    obs = [homozygous_variant, heterozygous, homozygous_reference]
    exp = hwe_expected(*obs)
    
    # Debug information
    st.write(f"Total Individuals: {sum(obs)}")
    st.write(f"Allele frequencies: p = {exp[2]**0.5}, q = {exp[0]**0.5}")

    # Display observed and expected counts
    st.write("### Observed and Expected Genotype Counts")
    results_df = pd.DataFrame({
        'Genotype': ['Homozygous Variant', 'Heterozygous', 'Homozygous Reference'],
        'Observed': obs,
        'Expected': exp
    })
    
    st.table(results_df)
    
    # Plot observed and expected counts
    st.write("### Observed vs. Expected Genotype Counts")
    fig, ax = plt.subplots()
    bar_width = 0.35
    index = np.arange(3)
    
    bar1 = ax.bar(index, obs, bar_width, label='Observed')
    bar2 = ax.bar(index + bar_width, exp, bar_width, label='Expected')
    
    ax.set_xlabel('Genotype')
    ax.set_ylabel('Count')
    ax.set_title('Observed vs. Expected Genotype Counts')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(['Homozygous Variant', 'Heterozygous', 'Homozygous Reference'])
    ax.legend()
    
    st.pyplot(fig)
    
    try:
        # Perform HWE chi-square test
        chi2, p_value = chisquare(obs, exp)
        st.write(f"### HWE Chi-Square Test")
        st.write(f"Chi-Square Statistic: {chi2:.4f}")
        st.write(f"P-Value: {p_value:.4f}")
    except Exception as e:
        st.write(f"An error occurred during HWE chi-square test: {e}")
    
    try:
        # Perform Fisher's Exact Test for 2x3 table
        if st.button("Perform Fisher's Exact Test for 2x3 Table"):
            # Create the 2x3 contingency table
            table = np.array([obs, exp])
            fisher_test = Table(table).test_nominal_association()
            p_value_fisher = fisher_test.pvalue
            st.write(f"### Fisher's Exact Test for 2x3 Table")
            st.write(f"P-Value: {p_value_fisher:.4f}")
    except Exception as e:
        st.write(f"An error occurred during Fisher's Exact Test for 2x3 Table: {e}")

    try:
        # Perform Fisher's Exact Test for 2x3 table using Monte Carlo simulation
        if st.button("Perform Fisher's Exact Test for 2x3 Table using Monte Carlo simulation"):
            # Convert counts to integers
            obs = [int(x) for x in obs]
            exp = [int(x) for x in exp]
            table = np.array([obs, exp])
            p_value_fisher = monte_carlo_fishers_exact(table)
            st.write(f"### Fisher's Exact Test for 2x3 Table using Monte Carlo Simulation")
            st.write(f"P-Value: {p_value_fisher:.4f}")
    except Exception as e:
        st.write(f"An error occurred during Fisher's Exact Test for 2x3 Table: {e}")

    try:
        # Perform permutation test for 2x3 table
        if st.button("Perform Permutation Test for 2x3 Table"):
            table = np.array([obs, exp])
            p_value_perm = permutation_test(table)
            st.write(f"### Permutation Test for 2x3 Table")
            st.write(f"P-Value: {p_value_perm:.4f}")
    except Exception as e:
        st.write(f"An error occurred during Permutation Test for 2x3 Table: {e}")

if page == "Allele Frequency Evolution":
    st.header("Allele Frequency Evolution Simulation")
    
    # User inputs
    p0 = st.slider("Initial Allele Frequency (p0)", 0.0, 1.0, 0.5)
    selection_coef = st.slider("Selection Coefficient (s)", 0.0, 1.0, 0.1)
    migration_rate = st.slider("Migration Rate (m)", 0.0, 1.0, 0.01)
    population_size = st.number_input("Population Size (N)", 1, 100000, 1000)
    generations = st.number_input("Number of Generations", 1, 1000, 100)
    
    # Run simulation
    if st.button("Run Simulation"):
        freqs = simulate_allele_freq(p0, selection_coef, migration_rate, population_size, generations)
        
        # Plot results
        plt.figure(figsize=(10, 6))
        plt.plot(freqs, label='Allele Frequency')
        plt.xlabel("Generations")
        plt.ylabel("Allele Frequency")
        plt.title("Allele Frequency Evolution Over Generations")
        plt.legend()
        st.pyplot(plt)

if page == "Statistical Power Analysis":
    st.header("Statistical Power Analysis")
    
    # User inputs
    effect_size = st.number_input("Effect Size (w)", 0.0, 2.0, 0.5)
    nobs = st.number_input("Sample Size (N)", 1, 100000, 1000)
    alpha = st.slider("Significance Level (alpha)", 0.0, 1.0, 0.05)
    
    # Calculate power
    if st.button("Calculate Power"):
        power = calculate_power(effect_size, nobs, alpha)
        st.write(f"Statistical Power: {power:.4f}")
        
        # Provide recommendations
        if power < 0.8:
            st.write("The statistical power is less than 0.8. Consider increasing the sample size for more reliable results.")
        else:
            st.write("The statistical power is sufficient for reliable results.")

        # Visualize power analysis
        sample_sizes = np.arange(10, nobs * 2, 10)
        powers = [calculate_power(effect_size, size, alpha) for size in sample_sizes]

        plt.figure(figsize=(10, 6))
        plt.plot(sample_sizes, powers, label='Power Curve')
        plt.axhline(0.8, color='r', linestyle='--', label='0.8 Power Threshold')
        plt.xlabel("Sample Size")
        plt.ylabel("Statistical Power")
        plt.title("Power Analysis Curve")
        plt.legend()
        st.pyplot(plt)

if page == "Clinical Interpretation":
    st.header("Clinical Interpretation of HWE Deviations")
    
    st.write("""
    In clinical genetics, deviations from Hardy-Weinberg Equilibrium (HWE) in a cohort of patients with a specific condition can indicate several things:
    - **Population Stratification**: Differences in allele frequencies due to population substructure.
    - **Selection**: The condition may be associated with certain genotypes, leading to an overrepresentation of these genotypes in the cohort.
    - **Genotyping Errors**: Technical errors in genotyping can lead to deviations.
    """)
    
    st.write("To help interpret deviations from HWE, use the tools below to analyze your data:")
    
    # User inputs for genotype counts
    st.write("Enter the observed genotype counts:")
    homozygous_variant = st.number_input("Homozygous Variant", 0, 100000, 0, key="ci_hom_var")
    heterozygous = st.number_input("Heterozygous", 0, 100000, 0, key="ci_het")
    homozygous_reference = st.number_input("Homozygous Reference", 0, 100000, 0, key="ci_hom_ref")
    
    obs = [homozygous_variant, heterozygous, homozygous_reference]
    exp = hwe_expected(*obs)
    
    # Perform HWE chi-square test
    chi2, p_value = chisquare(obs, exp)
    st.write(f"HWE Chi-Square Test: chi2 = {chi2:.4f}, p-value = {p_value:.4f}")
    
    # Interpretation guidelines
    if p_value < 0.05:
        st.write("The p-value is less than 0.05, indicating a significant deviation from HWE.")
        st.write("This could suggest population stratification, selection for certain genotypes, or genotyping errors.")
        if homozygous_variant > exp[0]:
            st.write("The observed number of homozygous variants is higher than expected, which might suggest selection for this genotype in the cohort.")
        else:
            st.write("The observed number of homozygous variants is not higher than expected based on HWE.")
    else:
        st.write("The p-value is greater than 0.05, indicating no significant deviation from HWE.")
        st.write("The observed genotype frequencies are consistent with those expected under random mating conditions.")
    
    # Provide HWE cutoff
    st.write("### HWE Cutoff for Clinical Interpretation")
    cutoff_hom_var = exp[0] * 1.5  # Arbitrary cutoff value; adjust as needed
    st.write(f"If the observed number of homozygous variants is greater than {cutoff_hom_var:.2f}, it may suggest selection for this genotype in the cohort.")
    
    # Variant information
    variant = st.text_input("Enter Variant ID (e.g., rs12345) for Additional Information:")
    if variant:
        info = fetch_variant_info(variant)
        if info:
            st.write(f"Variant Information for {variant}:")
            st.write(info)
        else:
            st.write(f"Information for {variant} not found.")
