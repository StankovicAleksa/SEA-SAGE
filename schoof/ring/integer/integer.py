from ring import Ring,RingElement
class IntegerRing(Ring):
    def add(self,el_1,el_2):
        return Integer(el_1.val+el_2.val)
    def mul(self,el_1,el_2):
        return Integer(el_1.val*el_2.val)
    def mod(self,el_1,el_2):
        return Integer(el_1.val%el_2.val)
    def neg(self,el):
        return Integer(-el.val)
    def one():
        return Integer(1)
    def zero():
        return Integer(0)
class Integer(RingElement):
    def __init__(self,value):
        RingElement.__init__(self,IntegerRing())
        self.val=value

