# ------------------------------------------------------------
# Before running this plotter, please have all relevant root files you
# are running have in their names the particle type (i.e.  muon, pion, gamma,
# kaon0L, electron, etc.), as well as the particle gun energy
#
# Ex: single_e-_100gev_CLD_REC.root
#
#
# To run a group of files, follow the format:
#      python plot_tuples.py /path/to/PFONtuple/files/*.root --outdir /path/to/output/directory/
#
# To run with selections, add --select
#
# This script is meant as an example for how to make various plots.
# Users can modify this at their discretion.
# ------------------------------------------------------------

import ROOT
import math
import re
import random
import sys
import ctypes
import argparse
import os
from collections import defaultdict
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ------------------------------------------------------------
# Region definition (theta in degrees)
# ------------------------------------------------------------

REGIONS = {
    "barrel":     (50.0, 130.0),
    "transition": (40.0, 50.0),
    "endcap":     (15.0, 40.0),
}

# ------------------------------------------------------------
# Selection Functions
#
# Any selection criteria the user chooses can be
# defined in this section. Currently we base our
# selection on the CLD paper
# ------------------------------------------------------------




def selection(ev):
    return cut1(ev) and cut2(ev)

def cut1(ev):
    select = False
    if not ev.okAng:
        return select
    else:
        select = True
    return select

def cut2(ev):
    select = False
    if not ev.okPt:
        return select
    else:
        select = True
    return select




# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

#Get energy from filename
def parse_energy(fname):
    m = re.search(r'(\d+)\s*gev', fname.lower())
    if not m:
        raise RuntimeError(f"Could not extract energy from filename: {fname}")
    return float(m.group(1))

#Get costheta 
def theta_deg(px, py, pz):
    p = math.sqrt(px*px + py*py + pz*pz)
    if p == 0:
        return None
    costh = pz / p
    costh = max(-1.0, min(1.0, costh))
    return math.degrees(math.acos(costh))

#Parse region based on inputted theta
def which_region(theta):
    for name, (tmin, tmax) in REGIONS.items():
        if tmin <= theta < tmax:
            return name
    return None

#Get single particle name/id from filename
def parse_particle_name(fname: str) -> str:
    n = os.path.basename(fname).lower()
    if any(k in n for k in ['gamma', 'photon', 'gam']):
        return 'gamma'
    if any(k in n for k in ['e-', 'electron']):
        return 'electron'
    if any(k in n for k in ['kaon0L', 'k0l', 'kaon0long', 'k0long', 'k_0^l', 'kl']):
        return 'kaon0l'
    if any(k in n for k in ['pi-', 'pim', 'pion-', 'pi+', 'pip', 'pion+', 'pion']):
        return 'pion'
    if any(k in n for k in ['mu-', 'mum', 'muon-', 'mu+', 'mup', 'muon+']):
        return 'muon'
    return 'particle'

def parse_particle_pdg(fname: str) -> int:
    n = os.path.basename(fname).lower()
    if any(k in n for k in ['gamma', 'photon', 'gam']):
        return 22
    if any(k in n for k in ['e-', 'electron']):
        return 11
    if any(k in n for k in ['kaon0L', 'k0l', 'kaon0long', 'k0long', 'k_0^l', 'kl']):
        return 130
    if any(k in n for k in ['pi-', 'pim', 'pion-']):
        return -211
    if any(k in n for k in ['pi+', 'pip', 'pion+']):
        return 211 
    if any(k in n for k in ['mu-', 'mu', 'muon']):
        return 13 

    return 0

#PID category for 2D histograms
def pid_category(pdg):

    apdg = abs(int(pdg))

    if apdg == 13:
        return "mu"
    if apdg == 22:
        return "gamma"
    if apdg == 11:
        return "e"
    if apdg in (211, 321, 2212):
        return "CH"
    if apdg in (130, 2112, 310):
        return "NH"
    else:
        return "other"



# ------------------------------------------------------------
# Metric Functions
# ------------------------------------------------------------

def pid_efficiency_from_tuple(fname, tree_name, select):
    truth_pdg = parse_particle_pdg(fname)
    if truth_pdg == 0:
        print(f"[warn] cannot infer truth PID from {fname}")
        return None

    f = ROOT.TFile.Open(fname)
    t = f.Get(tree_name)

    correct = 0
    total = 0

    for ev in t:
        mc = ev.mc_pdg
        reco = ev.pfo_pdg

        if mc != truth_pdg:
            continue

        if select == True and selection(ev) == False:
                continue

        total += 1
        if reco == truth_pdg:
            correct += 1

    f.Close()

    if total == 0:
        return None

    eff = correct / total
    err = (eff * (1 - eff) / total) ** 0.5
    print("Eff: ", eff)
    return eff, err


def pid_efficiency_vs_energy(files, select):
    pts = []
    
    for fname in files:
        E = parse_energy(fname)
        treename = "PFONtuple"
        res = pid_efficiency_from_tuple(fname, treename, select)
        particle = parse_particle_name(fname)
        if res is None:
            eff = 0.0
            err = 0.0
        else:
            eff, err = res

        pts.append((E, 100.0*eff, 0.0, 100.0*err))

        print("Processing:", fname)
        print("  Energy:", E)
        print("  eff =", eff)

    if len(pts) == 0:
        print("WARNING: no points filled!")
        return ROOT.TGraphErrors()
    pts = sorted(pts)

    gr = ROOT.TGraphErrors(len(pts))
    for i, (x, y, ex, ey) in enumerate(pts):
        gr.SetPoint(i, x, y)
        gr.SetPointError(i, ex, ey)

    gr.SetTitle("PID Efficiency vs Energy;E [GeV];Efficiency [%]")
    gr.SetMarkerStyle(21)
    gr.SetMarkerColor(ROOT.kMagenta+1)
    gr.SetLineWidth(2)

    return gr, particle

def draw_pid_vs_energy(files, outdir, select):
    c2 = ROOT.TCanvas("c_pid", "", 800, 700)
    gr_pid, particle = pid_efficiency_vs_energy(files, select)
    gr_pid.SetLineColor(ROOT.kBlue+1)
    gr_pid.SetMarkerColor(ROOT.kBlue+1)
    gr_pid.SetMarkerStyle(21)
    gr_pid.SetMarkerSize(1.2)
    gr_pid.Draw("AP")
    out_png_all = os.path.join(outdir, f"pid_eff_vs_energy_{particle}.png")
    c2.SaveAs(out_png_all)


def pid_eff_vs_costheta(files, tree_name, select):
    eff_by_E = {}
    for fname in files:
        E = parse_energy(fname)
        if E is None:
            print(f"[warn] could not parse energy from {fname}, skipping")
            
        eff = ROOT.TEfficiency(
            f"eff_costh_{int(E)}GeV",
            f"PID efficiency vs cos#theta (E={int(E)} GeV);cos#theta;Efficiency", 40, -1.0, 1.0)
        ROOT.SetOwnership(eff, False)
        f = ROOT.TFile.Open(fname)
        t = f.Get(tree_name)
        if not t:
            raise RuntimeError(f"No PFONtuple in {fname}")
        n_den = 0
        n_num = 0
        target_pdg = parse_particle_pdg(fname)
        particle = parse_particle_name(fname)
        for ev in t:

            if ev.mc_pdg != target_pdg:
                continue
            if select == True and selection(ev) == False:
                continue

            px = ev.mc_Px
            py = ev.mc_Py
            pz = ev.mc_Pz
            p = (px*px + py*py + pz*pz) ** 0.5
            if p <= 0:
                continue
        
            costh = pz / p
            n_den += 1
        
            passed = ev.has_reco_pfo and ev.pfo_pdg == target_pdg
            if passed:
                n_num += 1
            
            eff.Fill(passed, costh)
                    
        f.Close()
        
        #Feldman Cousins
        eff.SetStatisticOption(ROOT.TEfficiency.kFFC)
        eff.SetConfidenceLevel(0.95)

        print(f"[info] E={E:g} GeV: num={n_num}, den={n_den}, pdg={target_pdg}")
    
        eff_by_E[E] = eff

    return eff_by_E, particle

def draw_pid_eff_vs_costheta(eff_by_E, particle, outdir, outfile):

    if not eff_by_E:
        print("[warn] No efficiencies found, nothing to draw")
        return
    c_all = ROOT.TCanvas("c_pid_eff_costh_full",
                         "PID efficiency vs cos(#theta) - All Energies", 900, 700)
    c_all.SetGrid()

    leg = ROOT.TLegend(0.55, 0.18, 0.88, 0.40)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    frame = c_all.DrawFrame(-1.0, 0.5, 1.0, 1.05)
    frame.SetTitle("PID Efficiency vs cos#theta;cos#theta;Efficiency[%]")

    markers = [
        ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp,
        ROOT.kOpenCircle, ROOT.kOpenSquare
    ]

    colors = [
        ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
        ROOT.kMagenta + 1, ROOT.kCyan + 1
    ]
    
    for i, E in enumerate(sorted(eff_by_E.keys())):
        eff = eff_by_E[E]

        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]

        eff.SetLineColor(color)
        eff.SetMarkerColor(color)
        eff.SetMarkerStyle(marker)
        eff.SetMarkerSize(1.2)

        eff.Draw("P SAME")

        leg.AddEntry(eff, f"E = {E:g} GeV", "p")

    leg.Draw()
    
    out_png_all = os.path.join(outdir, f"pid_efficiency_vs_costheta_{particle}.png")
    c_all.SaveAs(out_png_all)
    outfile.cd()

    for E, eff in eff_by_E.items():
        cname = f"c_pid_eff_costh_{int(E)}GeV"
        c = ROOT.TCanvas(cname, cname, 800, 600)
        c.SetGrid()

        eff.Draw("AP")
        eff.SetTitle(f"PID Efficiency vs cos(#theta) - {E} GeV;cos#theta;Efficiency[%]")
        out_png = os.path.join(outdir, f"pid_efficiency_vs_costheta_{particle}_{int(E)}GeV.png")
        c.SaveAs(out_png)

        eff.Write()


    c_all.Close()
        

def fill_pid_matrices(tree, categories, select):

    ncat = len(categories)

    # Create empty histograms
    h_counts = ROOT.TH2D(
        "h_pid_counts",
        "PID Confusion Matrix;MC Category;Reco Category;Counts",
        ncat, 0, ncat,
        ncat, 0, ncat
    )
    h_norm = ROOT.TH2D(
        "h_pid_norm",
        "PID Confusion Matrix (normalized);MC Category;Reco Category;Counts",
        ncat, 0, ncat,
        ncat, 0, ncat
    )

    # Set bin labels
    for i, cat in enumerate(categories):
        h_counts.GetXaxis().SetBinLabel(i+1, cat)
        h_counts.GetYaxis().SetBinLabel(i+1, cat)
        h_norm.GetXaxis().SetBinLabel(i+1, cat)
        h_norm.GetYaxis().SetBinLabel(i+1, cat)

    # Fill count matrix
    for event in tree:
        mc_cat = pid_category(event.mc_pdg)

        reco_cat = "other"
        reco_cat = pid_category(event.pfo_pdg)
        if select == True and selection(event) == False:
            continue
        passed = event.has_reco_pfo

        xbin = categories.index(mc_cat)
        ybin = categories.index(reco_cat)

        if passed:
            h_counts.Fill(xbin, ybin)
            if reco_cat == "other":
                print(f"Reco pdg: {event.pfo_pdg}_")

    for i in range(1, ncat+1):
        row_sum = sum(h_counts.GetBinContent(i, j) for j in range(1, ncat+1))
        if row_sum == 0:
            continue
        for j in range(1, ncat+1):
            val = h_counts.GetBinContent(i, j) / row_sum
            h_norm.SetBinContent(i, j, val)

    return h_counts, h_norm
def draw_pid_matrices(h_counts, h_norm, outdir, energy):
    os.makedirs(outdir, exist_ok=True)
    c = ROOT.TCanvas("pid_ct", f"PID Confusion Matrix - {energy} GeV", 800, 600)
    c.SetRightMargin(0.15)
    h_counts.SetStats(1)
    h_counts.GetZaxis().SetTitle("Counts")
    h_counts.SetTitle(f"PID Confusion Matrix - {energy} GeV")
    h_counts.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(outdir, f"pid_conf_{energy:.0f}GeV_counts.png"))

    h_norm.SetStats(1)
    h_norm.SetMaximum(1.0)
    h_norm.SetMinimum(0.0)
    h_norm.GetZaxis().SetTitle("Counts - Normalized")
    h_norm.SetTitle(f"PID Confusion Matrix - {energy} GeV")
    h_norm.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(outdir, f"pid_conf_{energy:.0f}GeV_norm.png"))


def make_pid_confusion_matrix(tree, outfile, outdir, energy, select):

    categories = ["mu", "gamma", "e", "NH", "CH", "other"]

    h_counts, h_norm = fill_pid_matrices(tree, categories, select)

    draw_pid_matrices(h_counts, h_norm, outdir, energy)


def fit_resolution(hist):
    if hist.GetEntries() < 30:
        return None

    f = ROOT.TF1("g", "gaus",
                 hist.GetMean() - 2*hist.GetRMS(),
                 hist.GetMean() + 2*hist.GetRMS())

    f.SetParameters(hist.GetMaximum(), hist.GetMean(), hist.GetRMS())

    if hist.Fit(f, "RQSN") != 0:
        return None

    sigma = abs(f.GetParameter(2))
    esigma = f.GetParError(2)

    return 100.0 * sigma, 100.0 * esigma

def save_energy_histograms_grouped(hdict, fitdict, energy, outdir):

    c = ROOT.TCanvas(f"c_{energy:.0f}GeV", "", 1200, 400)
    c.Divide(3, 1)

    region_order = ["barrel", "transition", "endcap"]

    for i, region in enumerate(region_order, start=1):
        c.cd(i)

        h = hdict.get(region)
        f = fitdict.get(region)

        if not h or h.GetEntries() == 0:
            continue

        h.SetLineWidth(2)
        h.SetLineColor(ROOT.kAzure-3)
        h.Draw("HIST")

        if f:
            f.SetLineColor(ROOT.kRed+1)
            f.SetLineWidth(2)
            f.Draw("SAME")

        h.GetXaxis().SetTitle("(E_{reco}/E_{true})-1")
        h.GetYaxis().SetTitle("Events")
        h.SetTitle(f"{energy:.0f} GeV")

        first_bin = h.FindFirstBinAbove(0)
        last_bin  = h.FindLastBinAbove(0)
        x_min = h.GetXaxis().GetBinLowEdge(first_bin)
        x_max = h.GetXaxis().GetBinUpEdge(last_bin)
        center = (x_max + x_min) / 2.0
        width = x_max - x_min
        new_half_width = (width * 1.1) / 2.0
        #h.GetYaxis().SetRangeUser(0, h.GetYaxis().GetXmax()*1.3)
        h.GetYaxis().SetRangeUser(0, 1.3*h.GetMaximum())
        h.GetXaxis().SetRangeUser(center - new_half_width, center + new_half_width)        
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextSize(0.05)
        label.DrawLatex(0.15, 0.85, region)

    c.SaveAs(os.path.join(outdir, f"energy_response_{energy:.0f}GeV.png"))


def fit_resolution_v2(hist):
    n = hist.GetEntries()
    if n < 20:
        print("[warn] No fit available - too few entries")
        return None

    rms = hist.GetRMS()
    if rms <= 0:
        print("[warn] No fit available - invalid RMS")
        return None

    f = ROOT.TF1(
        f"gaus_{hist.GetName()}",
        "gaus",
        hist.GetMean() - 2*rms,
        hist.GetMean() + 2*rms
    )

    f.SetParameters(hist.GetMaximum(), hist.GetMean(), rms)

    status = hist.Fit(f, "RQSN")

    if status == 0:
        sigma  = abs(f.GetParameter(2))
        esigma = f.GetParError(2)
    else:
        sigma  = rms
        esigma = rms / math.sqrt(2*n)
        f = None  # no valid Gaussian

    return {
        "sigma": 100.0 * sigma,
        "esigma": 100.0 * esigma,
        "fit": f
    }


def fit_gauss_twostep(h, default_range=(-0.3, 0.3)):
    """
    Two-step Gaussian fit on 'h' without any rebinning:
      1) Fit in a general window (default_range clipped to axis).
      2) Refit in +/- 3*sigma around the first fit's mean.
    Returns dict with:
      ok, h (clone), fn (final TF1), xmin, xmax, sigp, esigp, mu, sigma, stage
    """
    if (not h) or h.GetEntries() < 10:
        return {'ok': False}

    # Work on a detached clone (keep original binning, no Rebin)
    hfit = h.Clone(h.GetName() + '_fitclone')
    hfit.SetDirectory(0)

    # Axis limits and initial window
    ax = hfit.GetXaxis()
    xlo_tot, xhi_tot = ax.GetXmin(), ax.GetXmax()
    x1 = max(default_range[0], xlo_tot)
    x2 = min(default_range[1], xhi_tot)
    if x2 <= x1:
        x1, x2 = xlo_tot, xhi_tot

    # 1) General fit
    f1 = ROOT.TF1('g1_{}'.format(id(hfit)), 'gaus', x1, x2)
    amp0 = max(hfit.GetBinContent(hfit.GetMaximumBin()), 1.0)
    mean0 = hfit.GetMean()
    rms0  = hfit.GetRMS() or (x2 - x1) / 6.0
    f1.SetParameters(amp0, mean0, rms0)
    r1 = hfit.Fit(f1, 'RQSN')
    if int(r1) != 0:
        f1.SetRange(xlo_tot, xhi_tot)
        f1.SetParameters(amp0, 0.0 if abs(mean0) > 0.5 else mean0, rms0)
        if int(hfit.Fit(f1, 'RQSN')) != 0:
            return {'ok': False}
    mu1 = f1.GetParameter(1)
    sg1 = abs(f1.GetParameter(2)) or (x2 - x1) / 6.0

    # 2) Second fit in +/- 3*sigma from the first fit
    w2_lo = max(xlo_tot, mu1 - 3.0 * sg1)
    w2_hi = min(xhi_tot, mu1 + 3.0 * sg1)
    if w2_hi - w2_lo < max(ax.GetBinWidth(1) * 3.0, 1e-3):
        # Window collapsed; use the first fit as the result
        sigp = 100.0 * abs(f1.GetParameter(2))
        esigp = 100.0 * f1.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f1, 'xmin': x1, 'xmax': x2,
                'sigp': sigp, 'esigp': esigp, 'mu': mu1, 'sigma': sg1, 'stage': 'first'}

    f2 = ROOT.TF1('g2_{}'.format(id(hfit)), 'gaus', w2_lo, w2_hi)
    f2.SetParameters(f1.GetParameter(0), mu1, sg1)
    r2 = hfit.Fit(f2, 'RQSN')
    if int(r2) == 0:
        sigp  = 100.0 * abs(f2.GetParameter(2))
        esigp = 100.0 * f2.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f2, 'xmin': w2_lo, 'xmax': w2_hi,
                'sigp': sigp, 'esigp': esigp, 'mu': f2.GetParameter(1),
                'sigma': abs(f2.GetParameter(2)), 'stage': 'second'}
    else:
        # Fall back to first fit
        sigp  = 100.0 * abs(f1.GetParameter(2))
        esigp = 100.0 * f1.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f1, 'xmin': x1, 'xmax': x2,
                'sigp': sigp, 'esigp': esigp, 'mu': mu1, 'sigma': sg1, 'stage': 'first'}


    
def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input ROOT files (Must run PFONTuple)"
    )
    parser.add_argument(
        "--outdir",
        default="plots",
        help="Output directory for all plots and histograms"
    )
    parser.add_argument(
        "--select",
        action="store_true",
        help="Enables selection cuts in all plots"
    )
    

    args = parser.parse_args(argv)

    outdir = args.outdir
    files = args.inputs
    select = args.select
    if select == True:
        print("Running with selections")
    else:
        print("Running without selections")
    ensure_dir(outdir)

    
    hists = defaultdict(dict)
    outfile = ROOT.TFile(os.path.join(outdir, "output_hists.root"), "RECREATE")

    eff_by_E, particle = pid_eff_vs_costheta(files, "PFONtuple", select)
    draw_pid_eff_vs_costheta(eff_by_E, particle, outdir, outfile)
    draw_pid_vs_energy(files, outdir, select)

    for fname in files:
        energy = parse_energy(fname)
        
        f = ROOT.TFile.Open(fname)
        t = f.Get("PFONtuple")
        if not t:
            raise RuntimeError(f"No PFONtuple in {fname}")

        make_pid_confusion_matrix(t, outfile, outdir, energy, select)

        for r in REGIONS:
            if r not in hists[energy]:
                hists[energy][r] = ROOT.TH1D(
                    f"h_resE_{r}_{energy:g}",
                    ";(E_{reco}/E_{true}-1;Events",
                    200, -1.0, 1.0
                )
                hists[energy][r].SetDirectory(0)

        for ev in t:

            if select and not selection(ev):
                continue
            
            if not ev.has_reco_pfo:
                continue

            th = theta_deg(ev.mc_Px, ev.mc_Py, ev.mc_Pz)
            if th is None:
                continue

            region = which_region(th)
            if region is None:
                continue


            hists[energy][region].Fill(ev.resE)

        f.Close()
    graphs = {}
    for region in REGIONS:
        g = ROOT.TGraphErrors()
        g.SetName(f"g_{region}")
        graphs[region] = g

    for energy in sorted(hists.keys()):
        fits_this_energy = {}        
        outfile.mkdir(f"{int(energy)}GeV")
        outfile.cd(f"{int(energy)}GeV")
        for region, h in hists[energy].items():

            if h.GetEntries() == 0:
                print(f"[warn] cannot infer truth PID from {h}")
                continue

            res = fit_resolution_v2(h)
            if not res:
                continue
            g = graphs[region]
            n = g.GetN()
            g.SetPoint(n, energy, res["sigma"])
            g.SetPointError(n, 0.0, res["esigma"])
            fit = res["fit"]
            fits_this_energy[region] = fit

            h.SetName(f"h_res_{region}_{int(energy)}GeV")
            h.Write()

            if fit is not None:
                fit.SetName(f"fit_res_{region}_{int(energy)}GeV")
                fit.Write()

        outfile.cd()
        save_energy_histograms_grouped(hists[energy], fits_this_energy, energy, outdir)
            
    c = ROOT.TCanvas("c", "", 800, 600)
    c.SetGrid()

    
    colors = {
        "barrel": ROOT.kBlue+1,
        "transition": ROOT.kGreen+2,
        "endcap": ROOT.kMagenta+1,
    }
    markers = {
        "barrel": 20,
        "transition": 21,
        "endcap": 22
    }

    legend = {
        "barrel": "Barrel - #theta: 50#circ-130#circ",
        "transition": "Transition - #theta: 40#circ-50#circ",
        "endcap": "Endcap - #theta: 15#circ-40#circ"
    }
    
    first = True
    leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)

    for region, g in graphs.items():
        g.SetMarkerStyle(markers[region])
        g.SetMarkerColor(colors[region])
        g.SetLineColor(colors[region])

        if first:
            g.Draw("APE")
            g.GetXaxis().SetTitle("Energy [GeV]")
            g.GetYaxis().SetTitle("#sigma(E)/E [%]")
            g.SetTitle("Energy Resolution for all Detector Regions")
            g.GetYaxis().SetRangeUser(0, g.GetYaxis().GetXmax()*1.3)
            first = False
        else:
            g.Draw("P SAME")

        leg.AddEntry(g, legend[region], "lp")

    leg.Draw()
    c.SaveAs(os.path.join(outdir, "resolution_vs_energy_by_region.png"))

    outfile.Write()
    outfile.Close()
    print("Histograms written")
    
if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
