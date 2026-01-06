import transitions.transitions as transitions

if __name__ == "__main__":
    # El valor de esta variable puede se cambiada para correr el programa sin usar la CLI
    # Para más información, consulta el archivo README.md
    activate_cli = False

    if activate_cli:
        import cli.cli as cli
        cli.run_cli()
    else:
        transitions.get_transitions("LI-7-H", n_transitions=4)