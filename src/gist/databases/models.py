"""File containing the instructions to build the GIST schema using
the SQLalchemy ORM
"""

from __future__ import annotations
from typing import List

from sqlalchemy import (Float, ForeignKey, Integer, LargeBinary, String)
from sqlalchemy.orm import (DeclarativeBase, relationship,
                            Mapped, mapped_column)

# Base
class Base(DeclarativeBase):
    pass



# Entry Tables
# =============================================================================
class Protein(Base):
    __tablename__ = "protein"

    # Primary Key
    ID

    # Core information of a protein
    UniProt
    Gene
    Sequence
