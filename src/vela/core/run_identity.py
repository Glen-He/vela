"""跨阶段运行目录使用的稳定身份规则。"""

import re

from vela.core.errors import VelaError

RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$")


def validate_run_id(run_id: str) -> None:
    """拒绝不能安全用作单层目录名的运行身份。"""
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise VelaError(
            "run_id must start with an alphanumeric character and contain only "
            "letters, digits, dot, underscore, or hyphen"
        )
