import os.path
import difflib

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("Could not find required dependencies (numpy, pandas) for test script.")
    raise

nlines_header = 30
refpath = os.path.join(os.path.dirname(__file__), "references")


def run_and_capture(command):
    """Run a command and capture its output."""
    from subprocess import run, PIPE

    res = run(command, universal_newlines=True, stdout=PIPE, stderr=PIPE)
    assert res.returncode == 0, " ".join(command) + " failed with error:\n{}".format(
        res.stderr
    )
    return res


def read_HiggsBounds_results(filepath):
    """Read a HiggsBounds_result.dat file into a pandas DataFrame"""
    df = pd.read_csv(
        filepath,
        index_col="n",
        comment="#",
        sep=r"\s+",
        names=[
            "n",
            "Mh1",
            "Mh2",
            "Mh3",
            "Mhplus1",
            "HBresult",
            "chan",
            "obsratio",
            "ncomb",
            "additional1",
            "additional2",
        ],
    ).dropna(how="all")
    df.index = df.index.astype(int, copy=False)
    return df


def diff_pd(df1, df2):
    """Identify differences between two pandas DataFrames"""
    assert (df1.columns == df2.columns).all(), "DataFrame column names are different"
    if any(df1.dtypes != df2.dtypes):
        "Data Types are different, trying to convert"
        df2 = df2.astype(df1.dtypes)
    if df1.equals(df2):
        return None
    else:
        # need to account for np.nan != np.nan returning True
        diff_mask = (df1 != df2) & ~(df1.isnull() & df2.isnull())
        ne_stacked = diff_mask.stack()
        changed = ne_stacked[ne_stacked]
        changed.index.names = ["id", "col"]
        difference_locations = np.where(diff_mask)
        changed_from = df1.values[difference_locations]
        changed_to = df2.values[difference_locations]
        return pd.DataFrame(
            {"from": changed_from, "to": changed_to}, index=changed.index
        )


def test_HiggsBounds_effC():
    """Test HiggsBounds using effective coupling input."""
    workdir = os.path.join("tests", "mhmodplus_effC")
    run_and_capture(
        ["./HiggsBounds", "LandH", "effC", "3", "1", os.path.join(workdir, "mhmod+_")]
    )
    comp = read_HiggsBounds_results(
        os.path.join(refpath, "mhmod+_HiggsBounds_results.dat.effC.ref")
    )
    res = read_HiggsBounds_results(
        os.path.join(workdir, "mhmod+_HiggsBounds_results.dat")
    )
    diff = diff_pd(comp.drop("chan", axis=1), res.drop("chan", axis=1))
    assert diff is None, "\n" + str(diff)


def test_HiggsBounds_hadr():
    """Test HiggsBounds using effective coupling input."""
    workdir = os.path.join("tests", "mhmodplus_hadr")
    run_and_capture(
        ["./HiggsBounds", "LandH", "hadr", "3", "1", os.path.join(workdir, "mhmod+_")]
    )
    comp = read_HiggsBounds_results(
        os.path.join(refpath, "mhmod+_HiggsBounds_results.dat.hadr.ref")
    )
    res = read_HiggsBounds_results(
        os.path.join(workdir, "mhmod+_HiggsBounds_results.dat")
    )
    diff = diff_pd(comp.drop("chan", axis=1), res.drop("chan", axis=1))
    assert diff is None, "\n" + str(diff)


def test_HBeffC():
    """Test the HBeffC example program"""
    run_and_capture("./HBeffC")
    with open(os.path.join(refpath, "HBeffC-output.dat.ref"), "r") as ifile:
        comp = ifile.readlines()
    with open("HBeffC-output.dat") as ifile:
        res = ifile.readlines()
    diff = difflib.Differ().compare(comp, res)
    delta = "".join(x for x in diff if not x.startswith("  "))
    assert len(delta) == 0, "Files differ:\n{}".format(delta)


def test_HBwithchannelrates():
    """Test the HBwithchannelrates example program"""
    run_and_capture("./HBwithchannelrates")
    colnames = [
        "index",
        "Mh1",
        "Mh2",
        "Mh3",
        "tbeta",
        "HBresult0",
        "chan0",
        "obsratio0",
        "ncombined0",
        "llh0",
        "HBresult1",
        "chan1",
        "obsratio1",
        "ncombined1",
        "llh1",
    ]
    ref = pd.read_csv(
        os.path.join(refpath, "Mh125_HBwithchannelrates.dat.ref"),
        sep=r"\s+",
        index_col=0,
        names=colnames,
    )
    df = pd.read_csv(
        "Mh125_HBwithchannelrates.dat", sep=r"\s+", index_col=0, names=colnames
    )
    diff = diff_pd(
        ref.drop(["chan0", "chan1"], axis=1), df.drop(["chan0", "chan1"], axis=1)
    )
    assert diff is None, "\n" + str(diff)


def test_HBwithLHClikelihood():
    """Test the HBwithLHClikelihood example program"""
    run_and_capture("./HBwithLHClikelihood")
    colnames = [
        "index",
        "Mh1",
        "Mh2",
        "Mh3",
        "tbeta",
        "HBresult",
        "chan",
        "obsratio",
        "ncombined",
        "HBresult_all",
        "chan_all",
        "obsratio_all",
        "ncombined_all",
        "Hindex",
        "M_av",
        "nc",
        "cbin",
        "llh_CMS8",
        "llh_exp_CMS8",
        "llh_CMS13",
        "llh_exp_CMS13",
        "llh_ALTAS13",
        "llh_exp_ATLAS13",
        "llh_ATLAS20",
        "llh_exp_ATLAS20",
    ]
    ref = pd.read_csv(
        os.path.join(refpath, "Mh125_HBwithLHClikelihood.dat.ref"),
        sep=r"\s+",
        index_col=0,
        names=colnames,
    )
    df = pd.read_csv(
        "Mh125_HBwithLHClikelihood.dat", sep=r"\s+", index_col=0, names=colnames
    )
    diff = diff_pd(
        ref.drop(["chan", "chan_all"], axis=1), df.drop(["chan", "chan_all"], axis=1)
    )
    assert diff is None, "\n" + str(diff)


def test_HBwithLEPlikelihood():
    """Test the HBwithLEPlikelihood example program"""
    run_and_capture("./HBwithLEPlikelihood")
    colnames = [
        "mh",
        "HBres",
        "channel",
        "obsratio",
        "ncombined",
        "chisq_wo_uncertainty",
        "chisq_w_uncertainty",
        "channel_chisq",
    ]
    ref = (
        pd.read_csv(
            os.path.join(refpath, "HBchisq-output.dat.ref"),
            sep=r"\s+",
            names=colnames,
            comment="#",
        )
        .dropna(how="all")
        .reset_index(drop=True)
    )
    df = (
        pd.read_csv("HBchisq-output.dat", sep=r"\s+", names=colnames, comment="#")
        .dropna(how="all")
        .reset_index(drop=True)
    )
    diff = diff_pd(ref.drop("channel", axis=1), df.drop("channel", axis=1))
    assert diff is None, "\n" + str(diff)
