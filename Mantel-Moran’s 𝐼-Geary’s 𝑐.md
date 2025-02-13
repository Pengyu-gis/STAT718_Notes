# Mantel’s Test

## **1. Introduction**
Mantel (1967) first described this test and used it to examine spatial patterns in the occurrence of leukemia. The Mantel test examines the correlation between two distance matrices from the same samples.  The distance matrices can be based on actual data (species abundances, environmental measurements, spatial coordinates) or hypothesized data (e.g., dummy variables coding for treatments). 
### **Applications**
- **Spatial Analysis**: Checking spatial autocorrelation in geographic data.
- **Genetic Analysis**: Comparing genetic distance and geographic distance.
- **Ecological Research**: Investigating species distribution based on environmental factors.

### **Key Idea**
If two observations are spatially close, their attribute values should also be similar. Mantel’s Test quantifies this relationship by comparing two distance matrices.

---

## **2. Mathematical Notation**
Let’s define the key components used in Mantel’s Test:

### **Observations**
- $s_i$ and $s_j$ represent spatial locations (e.g., latitude/longitude).
- $Z(s_i)$ represents an observed attribute at location $s_i$.

### **Distance Matrices**
1. **Spatial Distance Matrix** $W_{ij}$:
   - $W_{ij}$ represents the spatial distance between two points.
   - A common choice is the **Euclidean distance**:
     $$W_{ij} = \|(s_i, t_i) - (s_j, t_j)\|$$

2. **Attribute Difference Matrix** $U_{ij}$:
   - Measures the difference in attribute values between two locations.
   - Some common definitions:
     $$U_{ij} = |Z(s_i) - Z(s_j)|$$
     or squared differences:
     $$U_{ij} = (Z(s_i) - Z(s_j))^2$$

---

## **3. Mantel’s Test Statistic**
Mantel’s test statistic is computed as:

$$M_2 = \sum_{i=1}^{n} \sum_{j=1}^{n} W_{ij} U_{ij}$$

where:
- $W_{ij}$ represents the spatial distance between locations $i$ and $j$.
- $U_{ij}$ represents the attribute difference.

This measures the correlation between spatial distance and attribute similarity.

---

## **4. Derivation of the Formula**
### **Step 1: Compute Distance Matrices**
- Compute $W_{ij}$ for all location pairs.
- Compute $U_{ij}$ based on attribute differences.

### **Step 2: Compute Mantel’s Statistic $M_2$**
- Multiply corresponding elements of $W_{ij}$ and $U_{ij}$.
- Sum over all pairs.

### **Step 3: Perform a Randomization Test**
To test the significance of $M_2$, we:
1. Shuffle attribute values randomly.
2. Recalculate $M_2$.
3. Repeat multiple times (typically 999 permutations).
4. Compare the observed $M_2$ with the random distribution.

### **Step 4: Compute the Z-Score**
To assess statistical significance, compute:

$$Z_{\text{obs}} = \frac{M_2 - E(M_2)}{\sqrt{\text{Var}(M_2)}}$$

where:
- $E(M_2)$ is the mean of permuted values.
- $\text{Var}(M_2)$ is the variance of permuted values.

---

## **5. Interpretation of Results**
- **If p-value $< 0.05$** → Reject the null hypothesis → There is significant spatial correlation.
- **If p-value $> 0.05$** → Fail to reject the null hypothesis → No significant correlation.

### **Z-score Interpretation**
- $Z > 1.96$ or $Z < -1.96$ → Significant spatial dependence at a 5% level.
- $Z \approx 0$ → No correlation.

---

## **6. Summary**
| **Step** | **Explanation** |
|----------|---------------|
| **1. Construct Distance Matrices** | Compute $W_{ij}$ (spatial) and $U_{ij}$ (attribute). |
| **2. Compute $M_2$** | Sum of element-wise product of matrices. |
| **3. Perform Randomization** | Shuffle attributes and recompute $M_2$ multiple times. |
| **4. Compute Z-score & p-value** | Compare observed $M_2$ with randomized distribution. |
| **5. Interpret Results** | Assess significance of spatial correlation. |

---

## **7. Final Notes**
- **Mantel’s Test** is useful for detecting spatial correlation but assumes a meaningful distance metric.
- **It is computationally expensive** for large datasets.
- **It does not handle non-linear relationships well**.

---

Reference:
1. https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/mantel-test/
2. 
