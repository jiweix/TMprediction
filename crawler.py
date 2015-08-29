################################################################
#
# Simple crawler
#
# 06/23/15, Jiwei
# Modified on 08/24/2015 for acquiring data from http://wardlab.cbs.umn.edu/yeast/
# Certain function and parameters only work for this website
# Data is writen to YeastTM.dat file
################################################################
import urllib2

# MaxCount to avoid infinite loop
MAXCOUNT = 5000
# Data was collected from this site
thisurl = "http://wardlab.cbs.umn.edu/yeast/"


def getPage(url):
    """ return html content from a web page

    :param url: web page url

    :return: string: all the html content from that page
    """
    try:
        content = urllib2.urlopen(url).read()
    except:
        content = None
    finally:
        return content


def getNextContent(page, startS, endS):
    """ Given two strings, find the content between these two strings in the current page

    :param page: the content of current working web page
    :param startS: content start point
    :param endS: content end point

    :return: string: the content between start point and end point
    """
    # return the next element with given start and end string from the current page content
    if page is None:
        return None, -1
    if page.find(startS) == -1:
        return None, -1
    startPoint = page.find(startS) + len(startS)
    endPoint = page.find(endS, startPoint+1)
    return page[startPoint:endPoint], endPoint


def fetchUrls(page, mainSite):
    """ get all the urls from a single page

    :param page: the working web page, updated in the while loop
    :param mainSite: use this string to check is the url is still in the main site

    :return: list: all the urls in the current page
    """
    urls = []
    while True:
        if page is None:
            break
        url, endPoint = getNextContent(page, 'HREF="', '"')
        page = page[endPoint:]
        if url is not None:
            if 'http' not in url:
                url = mainSite + url
            else:
                continue
            # Avoid run to other sites
            if url not in urls and mainSite in url:
                urls.append(url)
        else:
            break
    return urls


def fetchContent(page):
    """ extract the useful data from a page

    :param page: the working web page

    :return: list: all content in <pre></pre>, that is all the amino acids sequences
    """
    if page is None:
        return None
    contents = []
    while True:
        content, endPoint = getNextContent(page, '<pre>', '</pre>')       
        if content is not None:
            contents.append(content)
            page = page[endPoint:]
        else:
            break
    return contents


def crawlSite(siteUrl):
    """ get all the html content from http://wardlab.cbs.umn.edu/yeast/
        This function only work with this website

    :param siteUrl: the url

    :return: list: contents from all the pages
    """
    # Initialize
    toCrawl = fetchUrls(getPage(siteUrl+'left.htm'), siteUrl)
    crawled = []
    pageContent = []
    count = 0
    # Crawling while there are still urls in toCrawl list
    while toCrawl:
        site1 = toCrawl.pop()
        if site1 not in crawled:
            page = getPage(site1)
            pageContent.append(fetchContent(page))
            urlList = fetchUrls(page, siteUrl+'yprotdb/')      
            for url in urlList:
                if url not in toCrawl:
                    if url not in crawled:
                        toCrawl.append(url)
            crawled.append(site1)
            print site1
        count += 1
        if count > MAXCOUNT:
            break
    return pageContent 


def extract(item):
    """ get all the pieces of sequences, one sequence would be either transmenbrane domain or not

    :param item: a pieces of web page content

    :return: list: all sequences with label
    """
    fragments = []
    while True:
        frag, endPoint = getNextContent(item, '<', '<')
        if frag is not None:
            if 'blue' in frag:
                label = '1'
            else:
                label = '0'
            seq = frag[frag.find('>')+1:].replace('\n', '').strip()
            if len(seq) >= 20:
                fragments.append(label+ ' ' +seq)
            item = item[endPoint-1:]
        else:
            break
    return fragments


if __name__ == '__main__':
    siteToCrawl = thisurl
    fileToWrite = 'YeastTM.dat'
    with open(fileToWrite, 'w') as f:
        for items in crawlSite(siteToCrawl):
            for item in items:
                for strings in extract(item):
                    f.write("%s\n" % strings.encode('utf8'))