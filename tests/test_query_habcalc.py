import numpy as np
from habitex import ArchiveExplorer
import pytest

def test_query_habcalc():
    explorer = ArchiveExplorer(optimistic=True)
    gj876 = explorer.query_exo(hostname='GJ 876')
    assert gj876['pl_type'][2] == 'Water'
    assert gj876['in_hz_cons'][2]

test_query_habcalc()

def test_query_paper():
    explorer = ArchiveExplorer(optimistic=False)
    rosenthal_table = explorer.query_exo(table='ps',paper='Rosenthal et al. 2021')
    assert len(rosenthal_table) == 9
    assert min(rosenthal_table['dec'].values) == pytest.approx(-1.164,abs=1e-3)
    assert min(rosenthal_table['ra'].values) == pytest.approx(10.207,abs=1e-3)

test_query_paper()