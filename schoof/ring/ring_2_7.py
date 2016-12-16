from abc import abstractmethod,ABCMeta
class Ring:
    __metaclass__=ABCMeta
    class RingElement:
        def add(self,other):
            Ring.add(self,other)
    @abstractmethod
    def add(self,el_1,el_2):
        pass
    @abstractmethod
    def mul(self,el_1,el_2):
        pass 
    @abstractmethod
    def mod(self,el_1,el_2):
        pass 
    @abstractmethod
    def one():
        pass 
    @abstractmethod
    def zero():
        pass


