from abc import abstractmethod,ABC
class Field(ABC):
    @abstractmethod
    def add(self,el_1,el_2):
        pass
    @abstractmethod
    def mul(self,el_1,el_2):
        pass 
    def div(self,el_1,el_2):
        return self.mul(el_1,self.invert(el_2))
    @abstractmethod
    def neg(self,el_1):
        pass
    @abstractmethod
    def invert(self,el):
        pass
    def sub(self,el_1,el_2):
        return self.add(el_1,self.neg(el_2))
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
    @abstractmethod
    def one(self):
        pass 
    @abstractmethod
    def eq(self,el_1,el_2):
        pass
    @abstractmethod
    def zero(self):
        pass

class FieldElement:
    def __init__(self,field):
        self.field=field
    def __add__(self,other):
        return self.field.add(self,other)
    def __mul__(self,other):
        return self.field.mul(self,other)
    def __mod__(self,other):
        return self.field.mod(self,other)
    def __div__(self,other):
        return self.field.div(self,other)
    def __neg__(self):
        return self.field.neg(self)
    def __eq__(self,other):
        return self.field.eq(self,other)
    def __invert__(self):
        return self.field.invert(self)
    def __sub__(self,other):
        return self.field.sub(self,other)
    def __xor__(self,other):
        return self.field.power(self,other)
