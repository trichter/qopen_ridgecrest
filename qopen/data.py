from obspy import read

#file:///home/eule/data/datasets/qopen_ridgecrest/data/events/38443535/APL.CI.HNN..2019.185.182759.38443535.ms
DATA = '../data/events/{}/{}.{}.{}.{}.*.ms'

import logging
log = logging.getLogger('data')

def get_data(network, station, location, channel, starttime, endtime, event):
    evid = str(event.resource_id).split('/')[-1]
    fname = DATA.format(evid, station, network, channel, location)
    stream = read(fname, 'MSEED')
    stream.trim(starttime, endtime)
    #print(starttime, endtime)
    #print(stream)
    if not stream:
        log.debug('no data for %s.%s', network, station)
        print('no data for %s.%s' % (network, station))
    return stream
