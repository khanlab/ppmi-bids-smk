rule download_from_csv:
    input: 'url_csv/{prefix}.csv'
    output: directory('raw_zips/{prefix}')
    shell: 
        'mkdir -p {output} && for url in `cat {input}`; do wget $url --directory-prefix={output}; done '
        

