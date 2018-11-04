def sim_snow_for_one_water_year(tavg,precip):
    todays_snow = 0 # start on October 1.
    snow_list = []
    i=0
    for day in tavg: 
        if np.isnan(tavg[i]) == True: # Skip February 29 if it is present. The nan was inserted by the function that reads the model data OR by the bias correction code. See "consistency check option."
            i+=1
            continue
        todays_snow = melt_one_day(todays_snow,tavg[i])
        todays_snow = accum_one_day(todays_snow,tavg[i],precip[i])
        snow_list.append(todays_snow)
        i+=1
    snow_array = np.array(snow_list)
    return snow_array

def melt_one_day(start_swe,tavg): # Hock 2003. 
    global melt_factor , melt_thresh_temperature   # These constants are different for each location, determined by calibration to SNOTEL data.
    if tavg <= melt_thresh_temperature : return start_swe
    else:
        swe_delta = (tavg - melt_thresh_temperature) * melt_factor
        end_swe = start_swe - swe_delta
        if end_swe < 0 : end_swe = 0
        return end_swe
        
def accum_one_day(start_swe,tavg,precip): # Dingman et al. 2002., Lutz et al. 2010.
    global precip_fraction, melt_thresh_temperature #Precip fraction is always .1667.
    if tavg <= melt_thresh_temperature : rain = 0.0
    elif tavg > (melt_thresh_temperature + 6) : rain = 1.0
    else:
        rain = precip_fraction * (tavg - melt_thresh_temperature)
    snow_increment = (1.0 - rain) * precip
    if snow_increment < 0 : snow_increment = 0
    end_swe = start_swe + snow_increment
    return end_swe
    