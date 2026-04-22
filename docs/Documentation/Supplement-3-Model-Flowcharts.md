# Supplement 3 — Model Flowcharts

## Main Simulation Loop

```mermaid
flowchart TD
    A[Initialize] --> B[Read params.txt + config]
    B --> C[Build geometry + grid]
    C --> D[Load forcing time series]
    D --> E{Warmup complete?}
    
    E -- No --> F[Hyd → Transport only]
    F --> E
    
    E -- Yes --> G[Main loop: t = WARMUP → MAXT]
    
    subgraph "Each timestep (Δt = DELTI; bundled example = 360 s)"
        G --> H[Hyd: Saint-Venant solver]
        H --> I[bgboundary: interpolate BCs]
        I --> J[Transport: TVD + Crank-Nicolson]
        J --> K{Biogeochemistry enabled?}
        K -- Yes --> L[Biogeo: reactions + carbonate]
        K -- No --> M[Skip]
        L --> N[updateSPM: sediment dynamics]
        M --> N
        N --> O{Output interval?}
        O -- Yes --> P[Write OUT/ files]
        O -- No --> Q[Next timestep]
        P --> Q
    end
```

## Biogeochemical Reaction Network

```mermaid
flowchart TD
    subgraph "Primary Production"
        PP[μ = μ_max × f(T) × f(I) × min(fN, fP, fSi)]
        PP --> Phy1[Phy1: Diatoms]
        PP --> Phy2[Phy2: Non-diatoms]
    end
    
    subgraph "Nutrient Uptake"
        Phy1 --> NO3[NO₃ uptake]
        Phy1 --> NH4u[NH₄ uptake]
        Phy1 --> Si_u[Si uptake]
        Phy2 --> NO3
        Phy2 --> NH4u
        Phy1 --> PO4u[PO₄ uptake]
        Phy2 --> PO4u
    end
    
    subgraph "Losses & Recycling"
        Phy1 --> MORT[Mortality → TOC]
        Phy2 --> MORT
        TOC_pool[TOC] --> ADEGRAD[Aerobic degradation]
        ADEGRAD --> NH4r[NH₄ release]
        ADEGRAD --> PO4r[PO₄ release]
        ADEGRAD --> DIC_up[DIC increase]
        ADEGRAD --> O2_down[O₂ consumption]
    end
    
    subgraph "Nitrogen Transformations"
        NH4r --> NIT[Nitrification: NH₄ → NO₃]
        NIT --> O2_down
        NO3 --> DENIT[Denitrification: NO₃ → N₂]
    end
    
    subgraph "Gas Exchange"
        O2_ex[O₂ air-water exchange]
        CO2_ex[CO₂ air-water exchange]
    end
    
    subgraph "Carbonate System"
        DIC_up --> pH_solve[pH solver: DIC + AT → pH, pCO₂]
        CO2_ex --> DIC_up
    end
```

## Calibration Stages

```mermaid
flowchart LR
    S1[Stage 1: Hydrodynamics] --> S2[Stage 2: Sediment]
    S2 --> S3[Stage 3: Biogeochemistry]
    S3 --> S4[Stage 4: Carbonate]
    
    S1 --- P1[Chézy, Rs, D₀, K_VDB]
    S2 --- P2[Erosion, deposition rates]
    S3 --- P3[μ_max, k_nit, k_ox, ...]
    S4 --- P4[piston_velocity_scale]
```

## File I/O Structure

```mermaid
flowchart LR
    subgraph "INPUT/"
        params[params.txt]
        config[config_input.txt]
        geom[Geometry/]
        forcing[Forcing/]
        trib[Tributaries/]
        valid[Validation/]
    end
    
    subgraph "Model"
        CGEM[C_GEM_Daily.exe]
    end
    
    subgraph "OUT/"
        hydro[velocity.csv/.bin, waterDepth.csv/.bin, ...]
        wq[Phy1.csv/.bin, NO3.csv/.bin, O2.csv/.bin, ...]
        flux[Flux_Advection_*.csv/.bin, Flux_Dispersion_*.csv/.bin]
        rates[Reaction_*.csv/.bin, Diag_*.csv/.bin]
        geo_out[geometry.csv]
        post[figures/ and tables/ from Python tools]
    end
    
    params --> CGEM
    config --> CGEM
    geom --> CGEM
    forcing --> CGEM
    trib --> CGEM
    CGEM --> hydro
    CGEM --> wq
    CGEM --> flux
    CGEM --> rates
    CGEM --> geo_out
    CGEM --> post
```
