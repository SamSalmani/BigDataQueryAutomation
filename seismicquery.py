from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import kilometers2degrees
import os
from obspy.io.sac.sactrace import SACTrace

os.chdir(r'C:\Users\hsalmanitehrani\Desktop\Thesis\Code\Networks\BK')  
st = UTCDateTime("2010-01-01T00:00:00")  
et = UTCDateTime("2021-01-01T00:00:00") 

client = Client('NCEDC', timeout = 600) # for station
clientEQ = Client("USGS") # for event

com_listOUT = "HH?,HN?" # broadband & strong-motion
pre_tw_sec = 90  # 90 seconds before the origin time
tw_sec = 5 * 60  # 5 min duration

minmag = 3.0
maxmag = 5.0

output_unit = "ACC" # m/s/s

network_name = "NC"

# f1<f2<f3<f4. effective band is between f2 and f3
f1, f2, f3, f4 = 0.01, 0.05, 40, 50  # 0.01, 0.05, 40, 50 Hz

# get all stations in a specific network
inv = client.get_stations(network=network_name, station="*", channel="*",
                             starttime=st, endtime=et,
                             minlatitude=37.07906, maxlatitude=38.37746, minlongitude=-122.916, 
                             maxlongitude=-121.384,level="response")

 
for i in range(len(inv[0])):
    catalog = clientEQ.get_events(starttime=st, endtime=et, minmagnitude=minmag, maxmagnitude=maxmag, 
                              latitude=inv[0][i].latitude, longitude=inv[0][i].longitude, minradius=(kilometers2degrees(0)), 
                              maxradius=(kilometers2degrees(200)))
    
    stdir = os.getcwd() + "/" + "BK/" + inv[0][i].code
    
    if not os.path.exists(stdir):
                os.makedirs(stdir)  
    
    for j in range(len(catalog)):
        # EVENT INFORMATION
        event = catalog[j]
        origin = event.origins[0]
        origin_time = origin.time
        evlat = origin.latitude
        evlon = origin.longitude
        evdp_km = origin.depth / 1000
        evmag = event.magnitudes[0].mag

        evid = event.origins[0]['extra']['dataid']['value']
        
        starttime = UTCDateTime(origin_time - pre_tw_sec)
        endtime = UTCDateTime(origin_time + tw_sec)
        
        try:
            invOUT = client.get_stations(network=network_name, station=inv[0][i].code, channel="*", location = "*",
                             starttime=starttime, endtime=endtime,
                             latitude=evlat, longitude=evlon,
                             level="response")
        except Exception:
            continue
        
        
        try:
            channelOUT = client.get_stations(network=network_name, station=inv[0][i].code, channel="*",location = "*",
                                 starttime=starttime, endtime=endtime,
                                 latitude=evlat, longitude=evlon,
                                 level="channel")
        except Exception:
            continue
        
        
        try:
            stOUT = client.get_waveforms(network=network_name, station=inv[0][i].code, channel="*", 
                                     location= "*", starttime=starttime, endtime=endtime)
        except Exception:
            continue
        
        stOUT2 = stOUT.copy()
        
        # FOR THE AVAILABLE DATA, IT IS CONVERTED TO SAC FILE FORMAT, WHICH IS LIGHTER THAN mseed FORMAT
        for k in range(len(stOUT2)):

            try:
                stOUT2[k].remove_response(inventory=invOUT, pre_filt=(f1, f2, f3, f4), output=output_unit)
            except ValueError:
                continue

            sac = SACTrace.from_obspy_trace(stOUT2[k])
            sac.lcalda = True

            netOUT = stOUT2[k].stats.network
            comOUT = stOUT2[k].stats.channel
            staOUT = stOUT2[k].stats.station
            locOUT = stOUT2[k].stats.location
        
            stid = netOUT + "." + staOUT + "." + locOUT + "." + comOUT

            try:
                sta_coordinate = invOUT.get_coordinates(stid, starttime)
            except Exception:
                continue
            
            sac.stla = sta_coordinate['latitude']
            sac.stlo = sta_coordinate['longitude']
            sac.stel = sta_coordinate['elevation']
            sac.evla = evlat
            sac.evlo = evlon
            sac.evdp = evdp_km

            sac.kuser0 = "ncss"
            sac.kevnm = str(evid)
            sac.mag = evmag
        
            evdir = stdir + "/" + (str)(evid) + "_M" + (str)(evmag)
            
            if not os.path.exists(evdir):
                os.makedirs(evdir)
            
            outsacfi = evdir + "/" + comOUT +'.sac'

            sac.write(outsacfi)
