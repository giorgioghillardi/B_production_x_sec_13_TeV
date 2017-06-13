create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins pt --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins pt --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins pt --vector mass_window
combine_syst --measure ratio --channel 1 --bins pt

create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins pt --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins pt --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins pt --vector mass_window
combine_syst --measure ratio --channel 2 --bins pt

create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins pt --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins pt --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins pt --vector mass_window
combine_syst --measure ratio --channel 4 --bins pt


create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins y --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins y --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 1 --bins y --vector mass_window
combine_syst --measure ratio --channel 1 --bins y

create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins y --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins y --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 2 --bins y --vector mass_window
combine_syst --measure ratio --channel 2 --bins y

create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins y --vector signal_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins y --vector cb_pdf
create_yield_eff_syst_vector --measure ratio --calculate 1 --channel 4 --bins y --vector mass_window
combine_syst --measure ratio --channel 4 --bins y