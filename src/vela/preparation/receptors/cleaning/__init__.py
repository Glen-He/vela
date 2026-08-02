"""方法无关 receptor-only 基础结构制备入口。"""

from vela.preparation.receptors.cleaning.conformers import (
    resolve_alternate_conformation,
)
from vela.preparation.receptors.cleaning.workflow import prepare_receptors

__all__ = ["prepare_receptors", "resolve_alternate_conformation"]
