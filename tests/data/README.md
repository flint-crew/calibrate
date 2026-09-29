# Test data

These files are copied unchanged from Flint's test data
(https://github.com/flint-crew/flint, `flint/data/tests` and `flint/data/models`,
BSD 3-Clause licence), at commit `52159a6`.

| File | Description |
| --- | --- |
| `SB39400.RACS_0635-31.beam0.small.ms.zip` | A cut-down ASKAP beam: 36 antennas, 288 channels, 3 timesteps, WEIGHT and SIGMA only. Its visibilities are NaN and flagged and its UVWs are zero, so `msgen.py` fills them in. |
| `SB39433.B1934-638.beam0.calibrate.bin` | Bandpass solutions made by `calibrate` |
| `SB38969.B1934-638.beam35.aocalibrate.bin` | Bandpass solutions with many NaN solutions, a case that once broke Flint |
| `1934-638.calibrate.txt` | Flint's model of PKS B1934-638 |
