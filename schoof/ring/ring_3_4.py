from abc import abstractmethod,ABC
class Ring(ABC):
    @abstractmethod
    def add(self,el_1,el_2):
        pass
    @abstractmethod
    def mul(self,el_1,el_2):
        pass 
    @abstractmethod
    def neg(self,el_1):
        pass
    def sub(self,el_1,el_2):
        return el_1+self.neg(el_2)
    @abstractmethod
    def one():
        pass 
    @abstractmethod
    def zero():
        pass
    r"""Calculates el^exp """
    def power(self,el,exp):
    # I choose here to avoid recursion for saving stack space
        sol=self.one()
        tmp=el
        while  exp!=0 :
            if ( exp & 1 != 0 ):
                sol*=tmp
            tmp*=tmp
            exp >>= 1
        return sol

class RingElement:
    def __init__(self,ring):
        self.ring=ring
    def __add__(self,other):
        return self.ring.add(self,other)
    def __mul__(self,other):
        return self.ring.mul(self,other)
    def __neg__(self):
        return self.ring.neg(self)
    def __sub__(self,other):
        return self.ring.sub(self,other)
    def __xor__(self,other):
        return self.ring.power(self,other)
    def __eq__(self,other):
        return self.ring.equal(self,other)
