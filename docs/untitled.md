You are correct, and I apologize for the oversight. The absorption term should indeed be included in the final form of the integrated radiative transport equation (IRTE). Let's properly derive the IRTE with the appropriate absorption term included.

### Derivation of the IRTE

Starting from the steady-state radiative transfer equation:

$$
\mathbf{\Omega} \cdot \nabla L(\mathbf{r}, \mathbf{\Omega}) + \kappa(\mathbf{r}) L(\mathbf{r}, \mathbf{\Omega}) = \int_{4\pi} \kappa_s(\mathbf{r}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}, \mathbf{\Omega}') \, d\mathbf{\Omega}' + S(\mathbf{r}, \mathbf{\Omega}),
$$

we integrate this equation along the direction $\mathbf{\Omega}$ from some initial position $\mathbf{r}'$ to $\mathbf{r}$. Parameterize the path by $s$, such that $\mathbf{r} = \mathbf{r}' + s \mathbf{\Omega}$.

### Step 1: Rewrite the transport equation in terms of $s$

$$\frac{d L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega})}{ds} + \kappa(\mathbf{r}' + s \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) = \int_{4\pi} \kappa_s(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}') \, d\mathbf{\Omega}' + S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}),
$$

where $\frac{d}{ds} = \mathbf{\Omega} \cdot \nabla$.

### Step 2: Integrate both sides from $s = 0$ to $s = s_f$

$$
\int_0^{s_f} \left( \frac{d L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega})}{ds} + \kappa(\mathbf{r}' + s \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \right) ds = \int_0^{s_f} \left( \int_{4\pi} \kappa_s(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}') \, d\mathbf{\Omega}' + S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \right) ds.
$$

### Step 3: Integrate the directional derivative term

Using the Fundamental Theorem of Calculus:

$$
\int_0^{s_f} \frac{d L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega})}{ds} \, ds = L(\mathbf{r}, \mathbf{\Omega}) - L(\mathbf{r}', \mathbf{\Omega}).
$$

### Step 4: Integrate the absorption term

$$
\int_0^{s_f} \kappa(\mathbf{r}' + s \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds.
$$

### Step 5: Integrate the scattering and source terms

$$
\int_0^{s_f} \int_{4\pi} \kappa_s(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}') \, d\mathbf{\Omega}' \, ds,
$$

$$
\int_0^{s_f} S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds.
$$

Putting it all together, we get:

$$
L(\mathbf{r}, \mathbf{\Omega}) - L(\mathbf{r}', \mathbf{\Omega}) + \int_0^{s_f} \kappa(\mathbf{r}' + s \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds = \int_0^{s_f} \int_{4\pi} \kappa_s(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}') \, d\mathbf{\Omega}' \, ds + \int_0^{s_f} S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds.
$$

Rearranging to solve for $L(\mathbf{r}, \mathbf{\Omega})$:

$$
L(\mathbf{r}, \mathbf{\Omega}) = L(\mathbf{r}', \mathbf{\Omega}) - \int_0^{s_f} \kappa(\mathbf{r}' + s \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds + \int_0^{s_f} \int_{4\pi} \kappa_s(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}', \mathbf{\Omega}) L(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}') \, d\mathbf{\Omega}' \, ds + \int_0^{s_f} S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds.
$$

### Expressing in the form $L(\mathbf{r}, \mathbf{\Omega}) = K L(\mathbf{r}', \mathbf{\Omega}') + S(\mathbf{r}, \mathbf{\Omega})$

We define $K$ as the product of the collision term $C$ and the transport term $T$.

#### Transport Term $T$

The transport term represents the attenuation of radiance due to absorption along the path:

$$
T(\mathbf{r}', \mathbf{\Omega}' \rightarrow \mathbf{r}, \mathbf{\Omega}) = \exp\left(-\int_0^s \kappa(\mathbf{r}' + s' \mathbf{\Omega}) \, ds'\right).
$$

#### Collision Term $C$

The collision term represents the scattering interactions within the medium:

$$
C(\mathbf{r}', \mathbf{\Omega}' \rightarrow \mathbf{r}, \mathbf{\Omega}) = \kappa_s(\mathbf{r}', \mathbf{\Omega}' \rightarrow \mathbf{r}, \mathbf{\Omega}).
$$

#### Kernel $K$

The kernel $K$ is the product of $C$ and $T$:

$$
K(\mathbf{r}', \mathbf{\Omega}' \rightarrow \mathbf{r}, \mathbf{\Omega}) = \kappa_s(\mathbf{r}', \mathbf{\Omega}' \rightarrow \mathbf{r}, \mathbf{\Omega}) \exp\left(-\int_0^s \kappa(\mathbf{r}' + s' \mathbf{\Omega}) \, ds'\right).
$$

### Final Form of IRTE

The final form of the integrated radiative transport equation (IRTE) is:

$$
L(\mathbf{r}, \mathbf{\Omega}) = \int_{\text{path}} \kappa_s(\mathbf{r}', \mathbf{\Omega}', \mathbf{\Omega}) \exp\left(-\int_0^s \kappa(\mathbf{r}' + s' \mathbf{\Omega}) \, ds'\right) L(\mathbf{r}', \mathbf{\Omega}') \, d\mathbf{r}' \, d\mathbf{\Omega}' + \int_0^s S(\mathbf{r}' + s \mathbf{\Omega}, \mathbf{\Omega}) \, ds.
$$

This equation now properly includes both the scattering and absorption terms, ensuring that the absorption effects are correctly accounted for in the integrated form.