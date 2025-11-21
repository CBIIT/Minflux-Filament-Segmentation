# MINFLUX Fiber Analysis – Fiji/ImageJ Macro

This repository contains a **Fiji/ImageJ macro** that takes 3D MINFLUX localization data (`X, Y, Z` in microns) and:

- Builds a 2D KDE-like localization density map (X–Y).
- Enhances filamentous structures using **Tubeness**.
- Detects fibers using **Ridge Detection**.
- Creates a robust **16-bit label image** with one label per fiber.
- Assigns each original 3D localization to a **fiber_id**.
- Computes **per-fiber, point-based 3D statistics** in µm:
  - 3D center (`x, y, z`).
  - Length (µm) along the main in-plane axis.
  - 3D “width” (µm) estimated from radial spread of localizations around the axis.

All steps run inside Fiji/ImageJ and produce both **QC images** and **CSV files** for downstream analysis.

---

## Table of Contents

- [1. Overview of the Pipeline](#1-overview-of-the-pipeline)
- [2. Algorithms and Methods](#2-algorithms-and-methods)
  - [2.1 Input parsing](#21-input-parsing)
  - [2.2 KDE-based density map (X–Y)](#22-kde-based-density-map-xy)
  - [2.3 Tubeness filtering](#23-tubeness-filtering)
  - [2.4 Ridge detection](#24-ridge-detection)
  - [2.5 Robust 16-bit fiber labeling](#25-robust-16-bit-fiber-labeling)
  - [2.6 Point → fiber assignment](#26-point--fiber-assignment)
  - [2.7 Per-fiber 3D statistics](#27-per-fiber-3d-statistics)
- [3. Input Format](#3-input-format)
- [4. Requirements](#4-requirements)
- [5. Parameters and Tuning](#5-parameters-and-tuning)
- [6. How to Run in Fiji](#6-how-to-run-in-fiji)
- [7. Output Files](#7-output-files)
  - [7.1 Intermediate images](#71-intermediate-images)
  - [7.2 CSV files](#72-csv-files)
- [8. Interpreting Per-Fiber Metrics](#8-interpreting-per-fiber-metrics)
- [9. Troubleshooting](#9-troubleshooting)
- [10. License](#10-license)

---

## 1. Overview of the Pipeline

Given a single CSV file of localizations with columns `(X, Y, Z)` in microns, the macro performs:

1. **CSV parsing**
   - Robust detection of delimiter (comma, tab, whitespace).
   - Optional header detection.
   - Skips comments and non-numeric rows.

2. **KDE-like 2D density map (X–Y)**
   - Bins localizations into a 2D grid at `pxPerUm` pixels per micron.
   - Applies Gaussian blur with `gaussSigmaUm` in µm to approximate a kernel density estimate (KDE).

3. **Tubeness filtering**
   - Normalizes the blurred density map.
   - Converts to 8-bit.
   - Applies the **Tubeness** filter to enhance line-like structures.

4. **Ridge Detection**
   - Runs the **Ridge Detection** plugin on the tubeness image.
   - Extracts fiber-like ridges as ROIs and adds them to the ROI Manager.
   - Optionally saves the ridge properties table.

5. **16-bit labeling of fibers**
   - Builds a 16-bit label image where each pixel value is:
     - `0` for background.
     - `1, 2, ..., N` for fibers (one label per ridge ROI).
   - Labeling is done by rasterizing each ROI onto a temporary mask and writing a unique value into the label image.

6. **Point → fiber mapping**
   - For each original localization:
     - Maps `(X, Y)` to image coordinates.
     - Reads the label from the 16-bit label image.
     - If off-line, searches a neighborhood to find the nearest labeled pixel.

7. **Per-fiber 3D statistics**
   - Uses only the original `(X, Y, Z)` data for points assigned to each fiber.
   - Computes:
     - 3D center (mean X, mean Y, median Z).
     - Fiber length in µm along the PCA-derived main axis in X–Y.
     - 3D width in µm from average radial distances around the axis.

8. **Outputs**
   - A set of **TIFF images** documenting intermediate steps.
   - Two main **CSV files**:
     - Per-point assignments (`*_with_fiber_id.csv`).
     - Per-fiber metrics (`*_FiberStats_points_um.csv`).

All outputs are stored in:

> `<parent_dir>/<base_name>_FIJI_OUT/`

where `<base_name>` is the input CSV filename without extension.

---

## 2. Algorithms and Methods

### 2.1 Input parsing

The macro:

- Opens the CSV with `File.openDialog(...)`.
- Normalizes line endings to `\n`.
- Skips:
  - Empty lines.
  - Lines starting with `#`.

**Delimiter detection**:

- If the first non-comment, non-empty line contains:
  - `\t` → `delimMode = "tab"`.
  - `,`  → `delimMode = "comma"`.
  - otherwise → whitespace (`"ws"`), with multiple spaces collapsed.

**Header detection**:

- First non-comment line is tokenized.
- If the first two tokens are numeric → **no header**.
- Otherwise → **header present** and data start on the next line.

**Data extraction**:

- The first three numeric tokens of each valid row are interpreted as:
  - `X` (µm), `Y` (µm), `Z` (µm).
- Non-numeric or incomplete rows are ignored.
- Data are stored in arrays `X[]`, `Y[]`, `Z[]`.

### 2.2 KDE-based density map (X–Y)

The macro constructs a 2D KDE-like image over X–Y:

1. Compute min/max of X and Y:

   ```text
   minX, maxX, minY, maxY
   
2. Expand by a padding distance:


Expand by a padding distance:

    padUm  = 3 * gaussSigmaUm
    xminUm = minX - padUm
    xmaxUm = maxX + padUm
    yminUm = minY - padUm
    ymaxUm = maxY + padUm

Convert to pixel dimensions using `pxPerUm`:

    widthPx  = floor((xmaxUm - xminUm) * pxPerUm + 1.5)
    heightPx = floor((ymaxUm - yminUm) * pxPerUm + 1.5)

Initialize a 32-bit float image `KDE_counts_32f` with zeros.

For each point `(X[i], Y[i])`, map to pixel indices:

    ix = round((X[i] - xminUm) * pxPerUm)
    iy = heightPx - 1 - round((Y[i] - yminUm) * pxPerUm)   // Y flip to ImageJ coords

If `(ix, iy)` is inside bounds, increment pixel value by `1.0`.

Apply Gaussian blur:

    sigmaPx = gaussSigmaUm * pxPerUm
    run("Gaussian Blur...", "sigma=" + sigmaPx)

This yields a smoothed 2D density map approximating a Gaussian KDE.

---

### 2.3 Tubeness filtering

To enhance filamentous structures:

1. Normalize the blurred KDE:

       I_norm = (I - min) / (max - min)

2. Convert to 8-bit:

   - Multiply by 255.
   - Convert to 8-bit.
   - Save as `<base_name>_kde_gray8.tif`.

3. Re-open `<base_name>_kde_gray8.tif` and run Tubeness:

       run("Tubeness", "sigma=1.0000 use")

4. Normalize the resulting 32-bit image again to [0, 1], then:

   - Multiply by 255.
   - Convert to 8-bit.
   - Save as `<base_name>_tubeness8bit.tif`.

The tubeness image highlights line-like structures and suppresses noise.

---

### 2.4 Ridge detection

Ridge detection is applied on `Tubeness_8bit`:

- The macro uses a pre-defined parameter string:

      line_width=3
      high_contrast=230
      low_contrast=30
      estimate_width
      extend_line
      displayresults
      add_to_manager
      method_for_overlap_resolution=SLOPE
      sigma=1.37
      lower_threshold=1.19
      upper_threshold=10.03
      minimum_line_length=0
      maximum=0

- The Ridge Detection plugin:
  - Identifies ridge-like curves (fibers).
  - Outputs a table of ridge properties (optionally saved).
  - Adds ROIs representing the ridges to the ROI Manager.

If no ROIs are found, the macro stops with a message:

> `Ridge Detection found no fibers.`

---

### 2.5 Robust 16-bit fiber labeling

To map fibers into an image:

1. Get the number of ROIs:

       nRois0 = roiManager("count")

2. Create a 16-bit image `FiberLabels_16u` initialized to 0, with the same width and height as the tubeness image.

3. For each ROI index `ri = 0..nRois0-1`:

   - Create a temporary 8-bit black image `tmpMask`.
   - Select ROI `ri` on `tmpMask`.
   - Set `Line Width...` to `labelLineWidthPx`.
   - Draw the ROI in white (`255`).
   - Threshold `tmpMask` from 1 to 255 and create a selection.
   - Add that selection as a new ROI to the ROI Manager.
   - Select this new ROI on `FiberLabels_16u` and set all selected pixels to:

         value = ri + 1

   - Remove the temporary ROI and close `tmpMask`.

The final `FiberLabels_16u` is a label image where:

- `0` = background.
- `1..nRois0` = individual fibers (one label per original ridge ROI).

---

### 2.6 Point → fiber assignment

Each original localization is mapped back to a fiber label:

1. For each point `(X[i], Y[i])`, compute `(ix, iy)` using the same coordinate transform as in the KDE step.

2. Read the label from `FiberLabels_16u`:

       lab = getPixel(ix, iy)

3. If `lab > 0`:
   - Assign:

         assignID[i] = floor(lab)

4. If `lab == 0`, use a local neighborhood search:

   - For radius `r = 1..nearestSearchRadiusPx`:
     - Evaluate a square window `[ix-r..ix+r] × [iy-r..iy+r]` (clipped to image bounds).
     - If any pixel `lab2 > 0` is found:
       - Set:

             assignID[i] = floor(lab2)

       - Stop searching.

Points that never find a label keep `assignID[i] = 0` (unassigned).

A sanity check is performed:

- Count how many points have `assignID > 0`.
- Print:

      [INFO] Assigned points: <assignedCount> / <nPts>

- If `assignedCount == 0`, the macro exits with a message suggesting to check `pxPerUm` and ridge parameters.

---

### 2.7 Per-fiber 3D statistics

Let `maxF = nRois0` and consider fiber IDs `fid = 1..maxF`.

Only points with `1 ≤ fiber_id ≤ maxF` are included.

#### Pass 1 – per-fiber sums and covariance terms

For each point with fiber `fid`:

- Increment count:

      cnt[fid]++

- Accumulate sums:

      sumXf[fid] += X[i]
      sumYf[fid] += Y[i]
      sumZf[fid] += Z[i]
      sumXX[fid] += X[i]^2
      sumYY[fid] += Y[i]^2
      sumXY[fid] += X[i]*Y[i]

#### PCA in X–Y and fiber center

For each fiber with `cnt[fid] > 0`:

- Means:

      mx = sumXf[fid] / cnt[fid]
      my = sumYf[fid] / cnt[fid]

- Center in X–Y:

      cx[fid] = mx
      cy[fid] = my

- Variances and covariance:

      varx = sumXX[fid]/cnt[fid] - mx^2
      vary = sumYY[fid]/cnt[fid] - my^2
      cov  = sumXY[fid]/cnt[fid] - mx*my

- Principal axis angle:

      theta = 0.5 * atan2(2*cov, varx - vary)

- Unit vectors:
  - Along-axis direction `(ux, uy)`:

        ux[fid] = cos(theta)
        uy[fid] = sin(theta)

  - Perpendicular direction `(vx, vy)`:

        vx[fid] = -uy[fid]
        vy[fid] =  cos(theta)

#### Pass 2 – length and Z storage

For each point assigned to fiber `fid`:

- Displacement from center:

      dx = X[i] - cx[fid]
      dy = Y[i] - cy[fid]

- Projection along main axis:

      s = dx*ux[fid] + dy*uy[fid]

- Track min and max `s`:

      smin[fid] = min(smin[fid], s)
      smax[fid] = max(smax[fid], s)

- Store Z values in a contiguous array `Zall` for each fiber (via offsets).

For each fiber, sort per-fiber Z values and compute median Z:

- Median Z:

      cz[fid] = median of Z for fiber fid

#### Pass 3 – 3D radial distances and width

For each point assigned to fiber `fid`:

- Decompose displacement:

      dx = X[i] - cx[fid]
      dy = Y[i] - cy[fid]
      dxy_perp = |dx*vx[fid] + dy*vy[fid]|   // distance in X–Y perpendicular to axis
      dz       = Z[i] - cz[fid]             // distance in Z from median

- 3D radial distance to fiber axis:

      r3D = sqrt(dxy_perp^2 + dz^2)

- Accumulate:

      sumRad[fid] += r3D

For each fiber:

- Length in µm from along-axis extent:

      length_um_pts = smax[fid] - smin[fid]

- 3D width in µm, defined as twice the mean radial distance:

      width3D_um_pts = 2.0 * (sumRad[fid] / cnt[fid])

These values, together with center coordinates and point counts, are written into the per-fiber table `FiberStats_points_um`.

---

### 3. Input Format

The macro expects a text/CSV file with at least 3 numeric columns:

- Columns 1–3 are interpreted as:
  - `X` – X coordinate in microns.
  - `Y` – Y coordinate in microns.
  - `Z` – Z coordinate in microns.

Additional columns (if present) are ignored.

Accepted features:

- Delimiters:
  - Comma-separated.
  - Tab-separated.
  - Whitespace-separated (multiple spaces tolerated).
- Header:
  - Optional.
  - Automatically detected by examining the first non-comment line.
- Comments:
  - Lines starting with `#` are ignored.
- Blank lines:
  - Ignored.

---

### 4. Requirements

- Fiji or ImageJ with macro support (ImageJ 1.53f or newer).
- Plugins:
  - Tubeness filter.
  - Ridge Detection plugin.
- Standard ImageJ commands:
  - Gaussian Blur, Threshold, Create Selection, Set..., etc.

---

### 5. Parameters and Tuning

Key parameters in the macro:

- `pxPerUm`  
  Pixels per micron (e.g. 10 → 0.1 µm per pixel). Controls the discretization of the KDE.

- `gaussSigmaUm`  
  Gaussian sigma in microns used for KDE smoothing. Too small → noisy; too large → fibers blur together.

- `labelLineWidthPx`  
  Stroke width (in pixels) when rasterizing ROIs into the label image. Larger values help ensure points land on labels.

- `nearestSearchRadiusPx`  
  Neighborhood radius (in pixels) for assigning points that fall just off the line.

- `ridgeParams`  
  String controlling Ridge Detection behavior (line width, contrast thresholds, etc.). Tune these if fibers are not detected or fragmented.

---

### 6. How to Run in Fiji

1. Open Fiji/ImageJ.
2. Open the macro (e.g. `Plugins → New → Macro`) and paste the code.
3. Save it as `MINFLUX_FiberAnalysis.ijm` (optional but recommended).
4. Run the macro (`Run` button in the macro editor).
5. When prompted, select your `RawLocs` CSV file (`X, Y, Z` in microns).
6. The macro:
   - Runs KDE → Tubeness → Ridge Detection → labeling → stats.
   - Writes outputs into:

         <parent_dir>/<base_name>_FIJI_OUT/

   - Prints log messages including assigned point counts.

---

### 7. Output Files

All outputs are written into:

    <parent_dir>/<base_name>_FIJI_OUT/

where `<base_name>` is the input CSV filename without extension.

#### 7.1 Intermediate images

- `<base_name>_counts_32f.tif`  
  32-bit raw counts image (binned localizations before smoothing).

- `<base_name>_kde_blur32f_sigma<gaussSigmaUm>um_counts.tif`  
  32-bit blurred KDE-like image.

- `<base_name>_kde_gray8.tif`  
  8-bit normalized KDE image (0–255).

- `<base_name>_tubeness32f_sigma1um.tif`  
  32-bit tubeness response image.

- `<base_name>_tubeness8bit.tif`  
  8-bit tubeness image used for Ridge Detection.

- `<base_name>_fiber_labels16u.tif`  
  16-bit label image:
  - 0 = background.
  - 1..N = fiber labels.

- `<base_name>_fiber_labels16u_autocontrast.tif`  
  16-bit label image with contrast enhanced for quick visualization.

#### 7.2 CSV files

- `<base_name>_RidgeDetectionResults.csv` (optional)  
  Ridge Detection result table (if the Results window exists). Contains properties for each detected ridge.

- `<base_name>_with_fiber_id.csv`  
  Per-point assignments. Columns:
  - `X` (µm)
  - `Y` (µm)
  - `Z` (µm)
  - `fiber_id`:
    - 0 = unassigned.
    - 1..N = fiber label.

- `<base_name>_FiberStats_points_um.csv`  
  Per-fiber statistics derived from original 3D points. Columns:

      fiber_id,
      n_points,
      center_x_um,
      center_y_um,
      center_z_um,
      length_um_pts,
      width3D_um_pts

  Interpretation:
  - `fiber_id`  
    Integer label matching `FiberLabels_16u`.

  - `n_points`  
    Number of localizations assigned to that fiber.

  - `center_x_um`, `center_y_um`  
    Mean X and Y coordinates (µm) of the fiber’s points.

  - `center_z_um`  
    Median Z (µm) of the fiber’s points.

  - `length_um_pts`  
    Fiber length in µm, computed from the extent of projections along the PCA-derived main axis in X–Y.

  - `width3D_um_pts`  
    3D fiber width (µm), defined as:

        2 × (mean 3D radial distance of points to the fiber axis)

    where the axis is defined by the fiber center and principal direction in X–Y, with radial distance measured in full 3D.

---

### 8. Interpreting Per-Fiber Metrics

Example row:

    fiber_id  n_points  center_x_um  center_y_um  center_z_um  length_um_pts  width3D_um_pts
    1         4727      -13.885      0.223        -0.011       1.704          0.325

- `fiber_id = 1`  
  First fiber in the label image.

- `n_points = 4727`  
  Strongly sampled fiber.

- `center_x_um`, `center_y_um`, `center_z_um`  
  3D position of the fiber center (in µm).

- `length_um_pts = 1.704`  
  Fiber length along its main in-plane axis.

- `width3D_um_pts = 0.325`  
  Approximate 3D thickness of the fiber, derived from the radial distribution of localizations.

---

### 9. Troubleshooting

**No fibers detected (Ridge Detection found no fibers)**

- Check the tubeness image `<base_name>_tubeness8bit.tif`.
- Adjust `ridgeParams` (e.g. thresholds, sigma).
- Tweak `pxPerUm` and `gaussSigmaUm` to change the scale of KDE and tubeness.

**No localizations assigned to fibers**

- Confirm `pxPerUm` matches the scale of your input coordinates (in µm).
- Inspect `<base_name>_fiber_labels16u_autocontrast.tif` and compare to the localization distribution.
- Increase `labelLineWidthPx` and/or `nearestSearchRadiusPx`.

**NaNs or strange values in downstream analysis**

- Verify your CSV parser:
  - Correct delimiter (`,`)?
  - Decimal separator (`.` vs `,`)?
- Inspect `*_FiberStats_points_um.csv` in a text editor to confirm numeric fields are present.

---

### 10. License

Add your preferred license text here, for example:

MIT License:

    MIT License

    Copyright (c) 2025 Adib Keikhosravi

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
