# SSLASSO 1.2.3 (2025-10-16)
- Faster updating rule has been added for when n > p
- Updates to RCpp functions (e.g. Calloc to R_Calloc)

# SSLASSO 1.2-2 (2019-12-13)

- Fixed: SSLASSO now uses `if(obj, "matrix")` instead of `class(obj) == "matrix"` in preparation for changes in R 4.0.0