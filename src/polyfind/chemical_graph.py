"""Covalent facts declared by the molecular producer, independent of geometry."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ChemicalBondGraph:
    atom_count: int
    edges: tuple

    def __post_init__(self):
        if type(self.atom_count) is not int or self.atom_count < 1:
            raise ValueError('positive integer chemical atom count required')
        normalized=[]; seen=set()
        for edge in self.edges:
            if len(edge)!=3:
                raise ValueError('chemical bond requires two atom indices and its order')
            first,second,order=edge
            if (type(first) is not int or type(second) is not int or first==second or
                    not 0<=first<self.atom_count or not 0<=second<self.atom_count):
                raise ValueError('distinct in-range integer chemical atom indices required')
            if (isinstance(order,bool) or not isinstance(order,(int,float)) or
                    not math.isfinite(order) or order not in (1,2,3)):
                raise ValueError('declared single, double or triple covalent order required')
            pair=tuple(sorted((first,second)))
            if pair in seen:
                raise ValueError('duplicate chemical bond')
            seen.add(pair); normalized.append((*pair,int(order)))
        object.__setattr__(self,'edges',tuple(normalized))

    @property
    def bonds(self):
        return [(first,second) for first,second,_ in self.edges]

    def identity(self):
        return dict(schema='explicit_ordered_covalent_graph_v1',atom_count=self.atom_count,
                    bond_orders_0idx=[list(edge) for edge in self.edges])
