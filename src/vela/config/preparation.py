"""受体审计和基础制备默认参数的严格组装。"""

from collections.abc import Mapping

from vela.config.values import assert_keys, boolean, number, string, table
from vela.preparation.receptors.models import (
    AltlocConfig,
    CrystalContactConfig,
    ReceptorAuditConfig,
    ReceptorPreparationConfig,
)


def parse_audit(source: Mapping[str, object]) -> ReceptorAuditConfig:
    """组装 `[audit.crystal_contacts]`。"""
    audit = table(source, "audit", path="")
    assert_keys(
        audit,
        allowed={"crystal_contacts"},
        required={"crystal_contacts"},
        path="audit",
    )
    contacts = table(audit, "crystal_contacts", path="audit")
    keys = {"distance_A", "min_occupancy", "include_hydrogens"}
    assert_keys(contacts, allowed=keys, required=keys, path="audit.crystal_contacts")
    return ReceptorAuditConfig(
        crystal_contacts=CrystalContactConfig(
            distance_A=number(contacts, "distance_A", path="audit.crystal_contacts"),
            min_occupancy=number(
                contacts, "min_occupancy", path="audit.crystal_contacts"
            ),
            include_hydrogens=boolean(
                contacts, "include_hydrogens", path="audit.crystal_contacts"
            ),
        )
    )


def parse_preparation(source: Mapping[str, object]) -> ReceptorPreparationConfig:
    """组装 `[preparation.altloc]`。"""
    preparation = table(source, "preparation", path="")
    assert_keys(
        preparation,
        allowed={"altloc"},
        required={"altloc"},
        path="preparation",
    )
    altloc = table(preparation, "altloc", path="preparation")
    keys = {"preferred_label"}
    assert_keys(altloc, allowed=keys, required=keys, path="preparation.altloc")
    return ReceptorPreparationConfig(
        altloc=AltlocConfig(
            preferred_label=string(
                altloc, "preferred_label", path="preparation.altloc"
            ),
        )
    )
