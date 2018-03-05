"""
Multi-ion container
"""
import plasmapy

import fiasco


class IonCollection(object):
    """
    Container for holding multiple ions. Instantiate with ions, elements, or another
    ion collection.

    Examples
    --------
    """
    
    def __init__(self, *args, **kwargs):
        self._ion_list = []
        for item in args:
            if isinstance(item, fiasco.Ion):
                self._ion_list.append(item)
            elif isinstance(item, fiasco.IonCollection):
                self._ion_list += item._ion_list
            else:
                raise TypeError(f'{item} has unrecognized type {type(item)}',
                                'and cannot be added to collection.')
        # TODO: check for duplicates
        assert all([all(self[0].temperature == ion.temperature) for ion in self]), (
            'Temperatures for all ions in collection must be the same.')
        
    def __getitem__(self, value):
        return self._ion_list[value]
    
    def __contains__(self, value):
        if type(value) is str:
            el, ion = value.split()
            if '+' in ion:
                ion = int(ion.strip('+')) + 1
            value = f'{plasmapy.atomic.atomic_symbol(el)} {ion}'
        elif isinstance(value, fiasco.Ion):
            value = value.ion_name
        return value in [i.ion_name for i in self._ion_list]
    
    def __add__(self, value):
        return IonCollection(self, value)
    
    def __radd__(self, value):
        return IonCollection(value, self)
