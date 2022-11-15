# Copyright 2022 Tom Eulenfeld, MIT license

from obspy import UTCDateTime as UTC
from obspy.core.event import Catalog, Event, Magnitude, Origin, ResourceIdentifier, Pick, Arrival, WaveformStreamID


def convert_events(fin, fout):
    with open(fin) as f:
        data = f.read()
    events = []
    stations = set()
    for line in data.splitlines():
        try:
            evid, date, time, mag, lon, lat, dep, relocated = line.split()
            otime = UTC(f'{date} {time}')
        except:
            continue
        picks = []
        arrivals = []
        try:
            with open(f'../data/phases/{evid}.phase') as f:
                data2 = f.read()
        except FileNotFoundError as ex:
            print(ex)
        else:
            for line in data2.splitlines():
                if evid in line:
                    _, _, _, datetime2, *_ = line.split()
                else:
                    try:
                        net, sta, cha, loc, _, _, _, phase, _, _, _, _, picktime = line.split()
                    except:
                        print(evid, 'break')
                        break
                    if loc == '--':
                        loc = ''
                    id_ = f'{net}.{sta}..{cha}'
                    ptime = UTC(datetime2) + float(picktime)
                    pick = Pick(time=ptime, waveform_id=WaveformStreamID(*id_.split('.')), phase_hint=phase)
                    arrival = Arrival(pick_id=pick.resource_id, phase=phase, time_weight=1)
                    picks.append(pick)
                    arrivals.append(arrival)
                    stations.add(f'{net}.{sta}')
        dep = float(dep) * 1000
        origin = Origin(time=otime, latitude=lat, longitude=lon, depth=dep, arrivals=arrivals)
        magnitude = Magnitude(mag=mag, magnitude_type='Ml')
        id_ = ResourceIdentifier(evid)
        event = Event(magnitudes=[magnitude], origins=[origin], picks=picks, resource_id=id_)
        events.append(event)
    events = Catalog(events)
    events.write(fout, 'CSZ', compression=True)
    return events, stations


# ev, stas = convert_events('../data/select55cat.txt', '../data/cat55.csz')
# print('stations', stas)

ev, stas = convert_events('../data/select55+1cat.txt', '../data/cat55+1.csz')
print('stations', stas)
