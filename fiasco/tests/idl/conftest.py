import pytest
# Do not require hissw for tests
try:
    import hissw
except ImportError:
    pass


@pytest.fixture(scope="session")
def idl_environment():
    if idl_available():
        return hissw.Environment(ssw_packages=["sdo/aia"], ssw_paths=["aia"])
    else:
        pytest.skip(
            "A working IDL installation is not available. You will not be able to run portions of the test suite."
        )


def idl_available():
    try:
        import hissw
        _ = hissw.Environment().run("")
        return True
    except Exception:
        return False
