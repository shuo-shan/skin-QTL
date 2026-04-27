# Parameters
n_qtls <- 97  # Total number of QTLs
x_overlaps <- 2  # Number of QTLs overlapping FOS binding sites

# Assuming an extremely low probability of overlap under null hypothesis
# This is a placeholder; in a real analysis, this should be based on additional data or assumptions.
p_expected <- x_overlaps / n_qtls

# Perform the binomial test
result <- binom.test(x_overlaps, n_qtls, p = p_expected, alternative = "greater")

# Print the result
print(result)


# Define the contingency table
contingency_table <- matrix(c(5, 0, 92, 97), nrow = 2, byrow = TRUE,
                            dimnames = list(c("QTLs", "Controls"),
                                            c("Overlap", "No Overlap")))

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)

# Print the result
print(fisher_test_result)
