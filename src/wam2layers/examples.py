from typing import Callable


def download_volta() -> None:
    print("downloading volta")  # TODO: implement


def download_eiffel() -> None:
    print("downloading eiffel")  # TODO: implement


AVAILABLE_CASES: dict[str, Callable[[], None]] = {
    "example-volta": download_volta,
    "example-eiffel": download_eiffel,
}
