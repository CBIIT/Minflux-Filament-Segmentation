// Fiji / ImageJ macro
// MINFLUX RawLocs -> KDE (10 px/µm) -> normalized -> 8-bit KDE (saved)
// -> Re-open 8-bit KDE -> Tubeness -> Normalize -> 8-bit
// -> Ridge Detection -> robust 16-bit labeling via per-ROI mask (no Fill on original ROI)
// -> Map points -> fiber_id
// -> POINT-BASED per-fiber stats in µm (length & 3D width from original points)

////////////////////// USER SETTINGS //////////////////////
inputPath = File.openDialog("Choose RawLocs CSV (X,Y,Z in microns)");
// inputPath = "/path/to/RawLocs.csv";

pxPerUm       = 10;          // 10 pixels per micron (pixel = 0.1 µm)
gaussSigmaUm  = 0.05;        // KDE Gaussian sigma (µm)
padUm         = 3 * gaussSigmaUm;

labelLineWidthPx      = 3;   // stroke thickness for labeling lines (pixels)
nearestSearchRadiusPx = 6;   // search radius (px) when a point lands off-line

ridgeParams = "line_width=3 high_contrast=230 low_contrast=30 estimate_width " +
              "extend_line displayresults add_to_manager method_for_overlap_resolution=SLOPE " +
              "sigma=1.37 lower_threshold=1.19 upper_threshold=10.03 minimum_line_length=0 maximum=0";
////////////////////////////////////////////////////////////

requires("1.53f");
if (inputPath=="") exit("No file chosen.");
setBatchMode(true);

// -------- Helpers (no ternary operators) --------
function maxOf(a,b){ if(a>b) return a; else return b; }
function minOf(a,b){ if(a<b) return a; else return b; }

function lastSepIndex(p){
  i=-1; j=-1;
  k=indexOf(p,"/");
  while(k!=-1){ i=k; k=indexOf(p,"/",k+1); }
  k=indexOf(p,"\\");
  while(k!=-1){ j=k; k=indexOf(p,"\\",k+1); }
  if(i>j) return i; else return j;
}

function getParentDir(p){
  k=lastSepIndex(p);
  if(k==-1) return getDirectory("current");
  else return substring(p,0,k+1);
}

function baseNameNoExtFromPath(p){
  k=lastSepIndex(p);
  if(k==-1) name=p; else name=substring(p,k+1);
  d=-1; pos=indexOf(name,".");
  while(pos!=-1){ d=pos; pos=indexOf(name,".",pos+1); }
  if(d==-1) return name; else return substring(name,0,d);
}

function isNumber(s){
  v=parseFloat(s);
  if(v==v) return true; else return false; // NaN check
}

function splitSmart(line,mode){
  if(mode=="tab")   return split(line,"\t");
  if(mode=="comma") return split(line,",");
  t=replace(line,"\t"," ");
  do{oldt=t; t=replace(t,"  "," ");} while(t!=oldt);
  while(lengthOf(t)>0 && substring(t,0,1)==" ") t=substring(t,1);
  while(lengthOf(t)>0 && substring(t,lengthOf(t)-1,lengthOf(t))==" ") t=substring(t,0,lengthOf(t)-1);
  if(t=="") return newArray(0);
  return split(t," ");
}

// -------- Paths --------
parentDir = getParentDir(inputPath);
baseName  = baseNameNoExtFromPath(inputPath);
outDir    = parentDir + baseName + "_FIJI_OUT" + File.separator;
File.makeDirectory(outDir);

// -------- Read CSV --------
whole = File.openAsString(inputPath);
whole = replace(whole,"\r\n","\n");
whole = replace(whole,"\r","\n");
lines = split(whole,"\n"); nlines = lines.length;

first=""; firstIdx=-1;
for(i=0;i<nlines;i++){
  t=trim(lines[i]);
  if(t!="" && substring(t,0,1)!="#"){ first=t; firstIdx=i; break; }
}
if(first=="") exit("No non-empty data line.");

delimMode="ws";
if(indexOf(first,"\t")!=-1) delimMode="tab";
else if(indexOf(first,",")!=-1) delimMode="comma";

tokens=splitSmart(first,delimMode);
if(tokens.length<2) exit("Bad first line.");

hasHeader=true;
if(isNumber(tokens[0]) && isNumber(tokens[1])) hasHeader=false;
if(hasHeader) startRow=firstIdx+1; else startRow=firstIdx;

count=0;
for(r=startRow;r<nlines;r++){
  t=trim(lines[r]); if(t==""||substring(t,0,1)=="#") continue;
  tok=splitSmart(t,delimMode); if(tok.length<3) continue;
  if(isNumber(tok[0])&&isNumber(tok[1])&&isNumber(tok[2])) count++;
}
if(count==0) exit("No numeric rows.");

X=newArray(count); Y=newArray(count); Z=newArray(count);
idx=0;
for(r=startRow;r<nlines;r++){
  t=trim(lines[r]); if(t==""||substring(t,0,1)=="#") continue;
  tok=splitSmart(t,delimMode); if(tok.length<3) continue;
  if(!(isNumber(tok[0])&&isNumber(tok[1])&&isNumber(tok[2]))) continue;
  X[idx]=parseFloat(tok[0]); Y[idx]=parseFloat(tok[1]); Z[idx]=parseFloat(tok[2]); idx++;
}
nPts=idx;
if(nPts==0) exit("No numeric data rows.");

// -------- KDE grid @ 10 px/µm --------
minX=X[0]; maxX=X[0]; minY=Y[0]; maxY=Y[0];
for(i=1;i<nPts;i++){
  if(X[i]<minX)minX=X[i];
  if(X[i]>maxX)maxX=X[i];
  if(Y[i]<minY)minY=Y[i];
  if(Y[i]>maxY)maxY=Y[i];
}
xminUm=minX-padUm; xmaxUm=maxX+padUm; yminUm=minY-padUm; ymaxUm=maxY+padUm;

widthPx  = floor((xmaxUm-xminUm)*pxPerUm + 1.5);
heightPx = floor((ymaxUm-yminUm)*pxPerUm + 1.5);
if(widthPx<10) widthPx=10;
if(heightPx<10) heightPx=10;

newImage("KDE_counts_32f","32-bit black",widthPx,heightPx,1);
countsID = getImageID();
for(i=0;i<nPts;i++){
  ix=round((X[i]-xminUm)*pxPerUm);
  iy=heightPx-1-round((Y[i]-yminUm)*pxPerUm);
  if(ix<0||iy<0||ix>=widthPx||iy>=heightPx) continue;
  setPixel(ix,iy, getPixel(ix,iy)+1.0);
}
selectImage(countsID);
saveAs("Tiff", outDir+baseName+"_counts_32f.tif");

// -------- Blur in place (provenance save) --------
selectImage(countsID);
sigmaPx = gaussSigmaUm * pxPerUm;
run("Gaussian Blur...", "sigma="+sigmaPx);
saveAs("Tiff", outDir+baseName+"_kde_blur32f_sigma"+d2s(gaussSigmaUm,4)+"um_counts.tif");

// -------- Normalize float copy --------
selectImage(countsID);
run("Duplicate...", "title=KDE_norm_32f");
normID = getImageID();
selectImage(normID);
getStatistics(a,mn,mi,ma,sd);
range=ma-mi; if(range<=0) range=1e-12;
run("Subtract...", "value="+mi);
run("Divide...",   "value="+range);

// -------- Save clean 8-bit KDE --------
selectImage(normID);
run("Duplicate...", "title=KDE_gray8_tmp");
grayID = getImageID();
selectImage(grayID);
run("Multiply...", "value=255");
setOption("ScaleConversions", false);
run("8-bit");
setOption("ScaleConversions", true);
rename("KDE_gray8");
saveAs("Tiff", outDir+baseName+"_kde_gray8.tif");

// -------- Re-open 8-bit KDE and run Tubeness (matches manual) --------
open(outDir + baseName + "_kde_gray8.tif");
tubSrcID = getImageID();
selectImage(tubSrcID);
run("Tubeness", "sigma=1.0000 use");
tub32ID = getImageID();
rename("Tubeness_32f");
saveAs("Tiff", outDir+baseName+"_tubeness32f_sigma1um.tif");

// -------- Normalize tubeness and make 8-bit --------
selectImage(tub32ID);
getStatistics(ta,tmn,tmi,tma,tsd);
tr = tma - tmi; if(tr<=0) tr=1e-12;
run("Subtract...", "value="+tmi);
run("Divide...",   "value="+tr);
run("Duplicate...", "title=Tubeness_8bit_tmp");
tub8ID = getImageID();
selectImage(tub8ID);
run("Multiply...", "value=255");
setOption("ScaleConversions", false);
run("8-bit");
setOption("ScaleConversions", true);
rename("Tubeness_8bit");
saveAs("Tiff", outDir+baseName+"_tubeness8bit.tif");

// -------- Ridge Detection (px units; provenance) --------
selectImage(tub8ID);
roiManager("Reset");
run("Ridge Detection", ridgeParams);
if (isOpen("Results")) saveAs("Results", outDir+baseName+"_RidgeDetectionResults.csv");

// -------- Robust labeling (16-bit) via per-ROI mask --------
roiManager("Show All");
nRois0 = roiManager("count");
if (nRois0==0) exit("Ridge Detection found no fibers.");

selectImage(tub8ID); getDimensions(w,h,c,z,t);
newImage("FiberLabels_16u", "16-bit black", w, h, 1);
labID = getImageID();
selectImage(labID);
run("Select None");

for (ri=0; ri<nRois0; ri++){
  // 1) Create temporary mask the same size
  newImage("tmpMask","8-bit black", w, h, 1);
  tmpID = getImageID();
  selectImage(tmpID);
  run("Select None");

  // 2) Apply ROI to mask and stroke it (works for all ROI types incl. composite)
  roiManager("Select", ri);
  run("Line Width...", "line="+labelLineWidthPx);
  setColor(255);
  run("Draw");

  // 3) Convert non-zero mask pixels to a selection
  setThreshold(1, 255);
  run("Create Selection");
  if (selectionType()!=-1) {
    // Send this selection via ROI Manager to the label image
    roiManager("Add");
    tmpIdx = roiManager("count") - 1;

    selectImage(labID);
    roiManager("Select", tmpIdx);
    run("Set...", "value="+(ri+1));
    run("Select None");

    // Remove the temporary ROI from the manager to keep indices stable
    roiManager("Select", tmpIdx);
    roiManager("Delete");
  }
  // 4) Close the mask
  selectImage(tmpID); close();
}

selectImage(labID);
saveAs("Tiff", outDir + baseName + "_fiber_labels16u.tif");

// -------- Map points -> fiber_id (with local search fallback) --------
assignID = newArray(nPts);
for(i=0;i<nPts;i++) assignID[i]=0;

selectImage(labID);
for(i=0;i<nPts;i++){
  ix=round((X[i]-xminUm)*pxPerUm);
  iy=heightPx-1-round((Y[i]-yminUm)*pxPerUm);
  if(ix<0||iy<0||ix>=widthPx||iy>=heightPx) continue;
  lab = getPixel(ix,iy);               // 16-bit label value (exact ID or 0)
  if (lab>0) assignID[i]=floor(lab);
  else {
    found=0;
    for(r=1; r<=nearestSearchRadiusPx && found==0; r++){
      x0=maxOf(ix-r,0); x1=minOf(ix+r,widthPx-1);
      y0=maxOf(iy-r,0); y1=minOf(iy+r,heightPx-1);
      for(yy=y0; yy<=y1 && found==0; yy++){
        for(xx=x0; xx<=x1; xx++){
          lab2=getPixel(xx,yy);
          if(lab2>0){ assignID[i]=floor(lab2); found=1; break; }
        }
      }
    }
  }
}

// -------- Save X,Y,Z with fiber_id --------
outDelim=",";
if(delimMode=="tab")   outDelim="\t";
if(delimMode=="comma") outDelim=",";

out1="X"+outDelim+"Y"+outDelim+"Z"+outDelim+"fiber_id\n";
for(i=0;i<nPts;i++)
  out1 += d2s(X[i],12)+outDelim+d2s(Y[i],12)+outDelim+d2s(Z[i],12)+outDelim+assignID[i]+"\n";
File.saveString(out1, outDir+baseName+"_with_fiber_id.csv");

// -------- POINT-BASED per-fiber stats in µm --------
maxF = nRois0; // valid IDs are 1..maxF

// Pass 1: counts and sums for covariance (per fiber)
cnt  = newArray(maxF+1);
sumXf= newArray(maxF+1); sumYf= newArray(maxF+1); sumZf= newArray(maxF+1);
sumXX= newArray(maxF+1); sumYY= newArray(maxF+1); sumXY= newArray(maxF+1);

for(i=0;i<nPts;i++){
  fid = assignID[i];
  if (fid<1 || fid>maxF) continue;
  xi=X[i]; yi=Y[i]; zi=Z[i];
  cnt[fid]++; sumXf[fid]+=xi; sumYf[fid]+=yi; sumZf[fid]+=zi;
  sumXX[fid]+=xi*xi; sumYY[fid]+=yi*yi; sumXY[fid]+=xi*yi;
}

// Means, PCA axes (per fiber)
cx = newArray(maxF+1); cy = newArray(maxF+1);
ux = newArray(maxF+1); uy = newArray(maxF+1);
vx = newArray(maxF+1); vy = newArray(maxF+1);

for(fid=1; fid<=maxF; fid++){
  n=cnt[fid]; if(n==0) continue;
  mx = sumXf[fid]/n; my = sumYf[fid]/n;
  cx[fid]=mx; cy[fid]=my;

  varx = sumXX[fid]/n - mx*mx;
  vary = sumYY[fid]/n - my*my;
  cov  = sumXY[fid]/n - mx*my;
  theta = 0.5 * atan2(2*cov, varx - vary);
  ux[fid] = cos(theta);  uy[fid] = sin(theta);
  vx[fid] = -uy[fid];    vy[fid] =  cos(theta);
}

// Pass 2: Z storage for median + along-axis span
offset = newArray(maxF+2);
offset[1] = 0;
for(fid=1; fid<=maxF; fid++) offset[fid+1] = offset[fid] + cnt[fid];
ZallLen = offset[maxF+1];
Zall = newArray(ZallLen);
filled = newArray(maxF+1);

smin = newArray(maxF+1); smax = newArray(maxF+1);
for(fid=1; fid<=maxF; fid++){ smin[fid]=1e12; smax[fid]=-1e12; }

for(i=0;i<nPts;i++){
  fid = assignID[i];
  if (fid<1 || fid>maxF) continue;
  dx = X[i] - cx[fid];
  dy = Y[i] - cy[fid];
  s  = dx*ux[fid] + dy*uy[fid];   // along-axis projection (µm)
  if (s<smin[fid]) smin[fid]=s;
  if (s>smax[fid]) smax[fid]=s;

  p = offset[fid] + filled[fid];
  Zall[p] = Z[i];
  filled[fid]++;
}

// Median Z per fiber
cz = newArray(maxF+1);
for(fid=1; fid<=maxF; fid++){
  n = cnt[fid]; if(n==0) { cz[fid]=0; continue; }
  start = offset[fid];
  Ztmp = newArray(n);
  for(k=0;k<n;k++) Ztmp[k] = Zall[start+k];
  Array.sort(Ztmp);
  if (n%2==1) cz[fid] = Ztmp[floor(n/2)];
  else        cz[fid] = 0.5 * (Ztmp[n/2 - 1] + Ztmp[n/2]);
}

// Pass 3: sum 3D radial distances to PCA axis (uses cz)
sumRad = newArray(maxF+1);
for(i=0;i<nPts;i++){
  fid = assignID[i];
  if (fid<1 || fid>maxF) continue;
  dx = X[i] - cx[fid];
  dy = Y[i] - cy[fid];
  dxy_perp = abs(dx*vx[fid] + dy*vy[fid]);   // µm
  dz = Z[i] - cz[fid];                       // µm
  sumRad[fid] += sqrt(dxy_perp*dxy_perp + dz*dz);
}

// ------- Emit point-based stats (Table + robust CSV write) -------
Table.create("FiberStats_points_um");
row=0;
csv2 = "fiber_id,n_points,center_x_um,center_y_um,center_z_um,length_um_pts,width3D_um_pts\n";

for(fid=1; fid<=maxF; fid++){
  n = cnt[fid]; if(n==0) continue;
  length_um_pts  = smax[fid] - smin[fid];       // µm
  width3D_um_pts = 2.0 * (sumRad[fid] / n);     // µm

  Table.set("fiber_id",       row, fid);
  Table.set("n_points",       row, n);
  Table.set("center_x_um",    row, cx[fid]);
  Table.set("center_y_um",    row, cy[fid]);
  Table.set("center_z_um",    row, cz[fid]);    // median Z
  Table.set("length_um_pts",  row, length_um_pts);
  Table.set("width3D_um_pts", row, width3D_um_pts);

  csv2 += fid + "," + n + "," +
          d2s(cx[fid],6) + "," + d2s(cy[fid],6) + "," + d2s(cz[fid],6) + "," +
          d2s(length_um_pts,6) + "," + d2s(width3D_um_pts,6) + "\n";
  row++;
}

// Robust write that does not depend on Table.save(name, path)
File.saveString(csv2, outDir + baseName + "_FiberStats_points_um.csv");

// Optional: quick preview of labels
selectImage(labID);
run("Enhance Contrast","saturated=0.35");
saveAs("Tiff", outDir+baseName+"_fiber_labels16u_autocontrast.tif");

setBatchMode(false);
print("[OK] Output saved to: "+outDir);
