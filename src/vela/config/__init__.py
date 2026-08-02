"""项目配置的公共入口。"""

from vela.config.loader import load_config
from vela.config.models import AppConfig, ConfigError, PathsConfig

__all__ = ["AppConfig", "ConfigError", "PathsConfig", "load_config"]
