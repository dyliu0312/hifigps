import pytest
from astropy import units as u
from hifigps.calculation import (
    dv2df, df2dv, dc_to_z, freq2z, z2freq, get_unit_factor,
    resolution, mainbeam, beam_solid_angle, drift_integration_time, 
    get_integration_time_FAST_HIIM, tb_sky, circular_area, SEFD, 
    flux_density_sensitivity, brightness_temperature_sensitivity,
    snu2tb, mass_to_brightness_temperature
)

def get_val(q):
    return q.value if hasattr(q, 'value') else q

def test_conversions():
    # dv2df
    df = dv2df(21 * u.km / u.s).to(u.MHz)
    assert pytest.approx(get_val(df), 0.01) == 0.099497
    
    # df2dv
    dv = df2dv(0.1 * u.MHz).to("km/s")
    assert pytest.approx(get_val(dv), 0.01) == 21.106
    
    # dc_to_z
    z = dc_to_z(100 * u.Mpc, 100 * u.Mpc, 0.1)
    assert pytest.approx(get_val(z), 0.01) == 0.11751
    
    # freq2z
    z_val = freq2z(1.29 * u.GHz)
    assert pytest.approx(get_val(z_val), 0.01) == 0.101089
    
    # z2freq
    freq = z2freq(0.1, 20 * u.km / u.s)
    assert pytest.approx(get_val(freq.to(u.Hz)), 100) == 1.29119965e+09
    
    # get_unit_factor
    assert get_unit_factor("cm", "m") == 0.01

def test_system_parameters():
    # resolution
    res = resolution()
    assert pytest.approx(get_val(res.to(u.arcmin)), 0.01) == 2.95067
    
    # mainbeam
    mb = mainbeam(0.1)
    assert pytest.approx(get_val(mb.to(u.arcmin)), 0.01) == 3.2457
    
    # beam_solid_angle
    bsa = beam_solid_angle(3 * u.arcmin).to(u.arcsec**2)
    assert pytest.approx(get_val(bsa), 0.1) == 36712.1
    
    # drift_integration_time
    dit = drift_integration_time()
    assert pytest.approx(get_val(dit.to(u.s)), 0.01) == 13.13
    
    # get_integration_time_FAST_HIIM
    it_hiim = get_integration_time_FAST_HIIM()
    assert pytest.approx(get_val(it_hiim.to(u.s)), 0.01) == 28.89

def test_sensitivity():
    # tb_sky
    sky = tb_sky()
    assert pytest.approx(get_val(sky.to(u.K)), 0.01) == 3.545
    
    # circular_area
    area = circular_area() * 0.7
    assert pytest.approx(get_val(area.to(u.m**2)), 1) == 49480
    
    # SEFD
    sefd = SEFD(20 * u.K).to(u.Jy)
    assert pytest.approx(get_val(sefd), 0.01) == 0.7812
    
    # flux_density_sensitivity
    fds = flux_density_sensitivity(t_sys=23.6 * u.K, delta_f=0.025 * u.MHz)
    assert pytest.approx(get_val(fds.to(u.mJy)), 0.001) == 0.8501

def test_mass_to_tb():
    # snu2tb
    z = 0.1
    t_rec = 20 * u.K
    freq = z2freq(z).to(u.GHz)
    t_sys = tb_sky(freq) + t_rec
    omega_mb = beam_solid_angle(mainbeam(z))
    sigma_sv = flux_density_sensitivity(t_sys)
    tb = snu2tb(sigma_sv, omega_mb, freq)
    assert pytest.approx(get_val(tb.to(u.mK)), 0.01) == 8.28
    
    # mass_to_brightness_temperature
    m2tb = mass_to_brightness_temperature(z, 1e9 * u.solMass).to(u.K)
    assert pytest.approx(get_val(m2tb), 0.01) == 10.795
