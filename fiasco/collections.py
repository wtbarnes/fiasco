"""
Multi-ion container
"""

import fiasco


class IonCollection(object):
    """
    Container for holding multiple ions. Can be instantiate with either a 
    """
    
    def __init__(self, *args, **kwargs):
        self._ion_list = []
        for item in args:
            if isinstance(item, fiasco.Ion):
                self._ion_list.append(item)
            elif isinstance(item, fiasco.Element):
                self._ion_list += [ion for ion in item]
            elif isinstance(item, type(self)):
                self._ion_list += item._ion_list
            else:
                raise TypeError('{} has an unrecognized type {} and cannot be added to collection.'
                                .format(item, type(item)))
        # TODO: check for duplicates
        # TODO: check all temperatures are the same
        
    def __getitem__(self, value):
        return self._ion_list[value]
    
    def __contains__(self, value):
        return value in [i.ion_name for i in self._ion_list]
    
    def __add__(self, value):
        return IonCollection(*(self._ion_list + [value]))
    
    def __radd__(self, value):
        return IonCollection(*([value] + self._ion_list))