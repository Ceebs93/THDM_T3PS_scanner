import os.path
import difflib
try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("Could not find required dependencies (numpy, pandas) for test script.")
    raise

refpath = os.path.join(os.path.dirname(__file__), 'references')


def run_and_capture(command):
    """Run a command and capture its output."""
    from subprocess import run, PIPE
    res = run(command, universal_newlines=True, stdout=PIPE, stderr=PIPE)
    assert res.returncode == 0, ' '.join(command)+" failed with error:\n{}".format(
        res.stderr)
    return res


def read_HiggsSignals_results(filepath):
    """Read a HiggsSignals_result.dat file into a pandas DataFrame"""
    df = pd.read_csv(filepath, index_col='n', comment='#', sep=r'\s+',
                     names=['n', 'Mh1', 'Mh2', 'Mh3', 'Mhplus1',
                            'csq_mu', 'csq_mh', 'csq_tot',
                            'nobs_mu', 'nobs_mh', 'nobs_tot',
                            'Pvalue', 'additional_1', 'additional_2']).dropna(how='all')
    df.index = df.index.astype(int, copy=False)
    return df


def diff_pd(df1, df2):
    """Identify differences between two pandas DataFrames"""
    assert (df1.columns == df2.columns).all(), \
        "DataFrame column names are different"
    if any(df1.dtypes != df2.dtypes):
        "Data Types are different, trying to convert"
        df2 = df2.astype(df1.dtypes)
    if df1.equals(df2):
        return None
    else:
        # need to account for np.nan != np.nan returning True
        diff_mask = ~np.isclose(df1, df2) & ~(df1.isnull() & df2.isnull())
        if not np.any(diff_mask):
            return None
        ne_stacked = diff_mask.stack()
        changed = ne_stacked[ne_stacked]
        changed.index.names = ['id', 'col']
        difference_locations = np.where(diff_mask)
        changed_from = df1.values[difference_locations]
        changed_to = df2.values[difference_locations]
        return pd.DataFrame({'from': changed_from, 'to': changed_to},
                            index=changed.index)


def test_HiggsSignals_effC():
    """Test HiggsSignals using effective coupling input."""
    workdir = os.path.join('tests', 'mhmodplus_effC')
    run_and_capture(['./HiggsSignals', 'latestresults', '2', 'effC', '3', '1',
                     os.path.join(workdir, 'mhmod+_')])
    comp = read_HiggsSignals_results(os.path.join(
        refpath, 'mhmod+_HiggsSignals_results.dat.effC.ref'))
    res = read_HiggsSignals_results(os.path.join(
        workdir, 'mhmod+_HiggsSignals_results.dat'))
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HiggsSignals_hadr():
    """Test HiggsSignals using hadronic input."""
    workdir = os.path.join('tests', 'mhmodplus_hadr')
    run_and_capture(['./HiggsSignals', 'latestresults', '2', 'hadr', '3', '1',
                     os.path.join(workdir, 'mhmod+_')])
    comp = read_HiggsSignals_results(os.path.join(
        refpath, 'mhmod+_HiggsSignals_results.dat.hadr.ref'))
    res = read_HiggsSignals_results(os.path.join(
        workdir, 'mhmod+_HiggsSignals_results.dat'))
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HS_2Higgses():
    """Test HS_2Higgses example program"""
    run_and_capture(['./HS_2Higgses'])

    def read_HS_2Higgses(file):
        return pd.read_csv(file, sep=r'\s+')

    comp = read_HS_2Higgses(os.path.join(refpath, '2Higgses_pdf2.dat.ref'))
    res = read_HS_2Higgses('results/2Higgses_pdf2.dat')
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HSeffC():
    """Test HSeffC example program"""
    run_and_capture(['./HSeffC'])

    def read_HSeffC(file):
        return pd.read_csv(file, comment='#', sep=r'\s+',
                           names=['mh', 'kappaF', 'kappaV',
                                  'Chis_mu', 'Chisq', 'ndf',
                                  'Htogaga_rate', 'HtoVV_rate', 'HtoFF_rate'])

    comp = read_HSeffC(os.path.join(refpath, 'HSeffC.dat.ref'))
    res = read_HSeffC('results/HSeffC.dat')
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HShadr():
    """Test HShadr example program"""
    run_and_capture(['./HShadr'])

    def read_HShadr(file):
        return pd.read_csv(file, comment='#', sep=r'\s+',
                           names=['mh', 'scale_ggf', 'scale_VH',
                                  'Chisq_mu', 'Chisq_mh', 'Chisq', 'ndf', 'PValue'])

    comp = read_HShadr(os.path.join(refpath, 'HShadr.dat.ref'))
    res = read_HShadr('results/HShadr.dat')
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HSwithSLHA():
    """Test HSwithSLHA example program"""
    run_and_capture(['./HSwithSLHA', '1', '../tests/SLHA_FHexample.fh'])

    def read_HSwithSLHA(file):
        return pd.read_csv(file, index_col=0, comment='#', sep=r'\s+',
                           names=['n', 'chisq_mu', 'chisq_mass', 'ndf'])

    comp = read_HSwithSLHA(os.path.join(
        refpath, 'SLHA_FHexample.fh-fromHS.ref'))
    res = read_HSwithSLHA('../tests/SLHA_FHexample.fh-fromHS')
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)


def test_HS_mass():
    """Test HS_mass example program"""
    run_and_capture(['./HS_mass'])

    def read_HS_mass(file):
        return pd.read_csv(file, comment='#', sep=r'\s+',
                           names=['mh', 'dmh', 'Chisq_mu', 'Chisq_mh', 'Chisq',
                                  'nassigned', 'ndf', 'Lambda'])

    for pdf in ['pdf1', 'pdf2', 'pdf3']:
        comp = read_HS_mass(os.path.join(refpath, 'HS_mass_'+pdf+'.dat.ref'))
        res = read_HS_mass('results/HS_mass_'+pdf+'.dat')
        diff = diff_pd(comp, res)
        assert diff is None, 'in: '+pdf+'\n'+str(diff)

def test_HS_SM_LHCRun1():
    """Test HS_SM_LHCRun1 example program"""
    run_and_capture(['./HS_SM_LHCRun1'])

    def read_HS_SM_LHCRun1(file):
        return pd.read_csv(file, comment='#', sep=r'\s+',
                           names=['mh', 'dmh', 'Chisq_mu', 'Chisq_mh', 'Chisq',
                                  'nassigned', 'ndf', 'Lambda'])

    for pdf in ['pdf1', 'pdf2', 'pdf3']:
        comp = read_HS_SM_LHCRun1(os.path.join(refpath, 'HS_SM_LHCrun1_mass_'+pdf+'.dat.ref'))
        res = read_HS_SM_LHCRun1('results/HS_SM_LHCrun1_mass_'+pdf+'.dat')
        diff = diff_pd(comp, res)
        assert diff is None, 'in: '+pdf+'\n'+str(diff)

def test_HSLimitEffC():
    """Test HSLimitEffC example program"""
    run_and_capture(['./HSLimitEffC'])

    def read_HSLimitEffC(file):
        return pd.read_csv(file, comment='#', sep=r'\s+')

    comp = read_HSLimitEffC(os.path.join(refpath, 'HSLimitEffC.dat.ref'))
    res = read_HSLimitEffC('results/HSLimitEffC.dat')
    diff = diff_pd(comp, res)
    assert diff is None, '\n'+str(diff)
