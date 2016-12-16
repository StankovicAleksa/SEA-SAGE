import sys
if ( sys.version_info >= (3,4) ):
    from ring.ring_3_4 import Ring,RingElement
elif ( sys.version_info >= (3,0) ):
    from ring.ring_3_0 import Ring,RingElement
elif ( sys.version_info >= (2,7) ):
    from ring.ring_2_7 import Ring,RingElement
else:
    raise Exception("You are using too old version of python. Please upgrade your python installation(>=2.7)")
